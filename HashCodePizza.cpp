#define RUN_ME_WITH_SH /*
g++ -Wall -std=c++14 -g -DIL_STD -DVERBOSE -o ${0%.*} ${0}          \
	-I/opt/ibm/ILOG/CPLEX_Studio1261/cplex/include                       \
	-I/opt/ibm/ILOG/CPLEX_Studio1261/concert/include                     \
	-L/opt/ibm/ILOG/CPLEX_Studio1261/cplex/lib/x86-64_linux/static_pic   \
	-L/opt/ibm/ILOG/CPLEX_Studio1261/concert/lib/x86-64_linux/static_pic \
	-lilocplex     \
	-lcplex        \
	-lcplexdistmip \
	-lconcert      \
	-lpthread      \
	-lm
exit
*/
#include <ilcplex/ilocplex.h>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <queue>

using namespace std;

int minNumCross = 3;
int maxRectSize = 12;

enum BlockType { RegularBlock = 0
               , HamBlock
               , SelectedRegularBlock
               , SelectedHamBlock
               , DoNotExistsBlock
};

char const * blockStrings[] {
	"-", "+", "x", "o", "#"
};

struct Base {
	Base() {   }
	Base( vector< BlockType > && b, int w, int h )
		: blocks{ move( b ) }, width{ w }, height{ h } {
	}

	void addBlock( BlockType type ) {
		switch ( type ) {
			case DoNotExistsBlock:     {   blocks.push_back( DoNotExistsBlock );       ++numDoNotExists;   break;   }
			case SelectedRegularBlock: {   blocks.push_back( SelectedRegularBlock );   ++numRegular;       break;   }
			case RegularBlock:         {   blocks.push_back( RegularBlock );           ++numRegular;       break;   }
			case SelectedHamBlock:     {   blocks.push_back( SelectedHamBlock );       ++numHam;           break;   }
			case HamBlock:             {   blocks.push_back( HamBlock );               ++numHam;           break;   }
		}
	}

	vector< BlockType > blocks;
	int width, height;
	int numRegular = 0;
	int numHam     = 0;
	int numDoNotExists = 0;

	int size() const {   return (int)blocks.size();   }
	int index( int i, int j ) const  {   return i + width*j;   }
	BlockType access( int i, int j ) const {   return blocks[index( i, j )];   }
};

Base LoadBaseFromFile( char const * fileName ) {
	FILE * fin = fopen( fileName, "rb" );
	Base base;
	if ( fscanf( fin, "%d %d %d %d\n", &base.height, &base.width, &minNumCross, &maxRectSize ) != 4 ) {
			cerr << "Invalid first line in input file" << endl;
			exit( 1 );
	}
	for ( int i = 0; i < (base.width+1)*base.height; ++i ) {
		char rc;
		if ( fscanf( fin, "%c", &rc ) != 1 ) {
			cerr << "Missing characters in input file" << endl;
			exit( 1 );
		}
		switch ( rc ) {
			case '\n': {                                         continue;    }
			case  '-': // fall through
			case  'T': {   base.addBlock( RegularBlock         );   break;    }
			case  '+': // fall through
			case  'H': {   base.addBlock( HamBlock             );   break;    }
			case  'x': {   base.addBlock( SelectedRegularBlock );   break;    }
			case  'o': {   base.addBlock( SelectedHamBlock     );   break;    }
			case  '#': {   base.addBlock( DoNotExistsBlock     );   break;    }
			default:   {
				cerr << "Unknown character: '" << rc << "' in input file" << endl;
				exit( 1 );
			}
		}
	}
	fclose( fin );
	return base;
}

class SubProblem {
	public:
		SubProblem( Base const & b )
			: env{}
			, model{ env }
			, cplex{ model }
			, base{ b } {

			InitBaseVariables();
			InitRectVariables();
			InitOverlapConstraints();
			InitObjectiveFunction();

			SetCplexParameters( cplex );
		}

		void Solve() {
			if ( !cplex.solve() ) {
				if ( cplex.getStatus() != IloAlgorithm::Infeasible ) {
					clog << "Optimization error, CPLEX status code: " << cplex.getStatus() << endl;
				} else {
					clog << "Infeasible, CPLEX status code: " << cplex.getStatus() << endl;
				}
			} else {
				FillSolvedBase();
				OutputSolution( cout );
			}
		}

		double GetBestScore() const {
			return cplex.getBestObjValue();
		}

		void PrintConstraints( ostream & out = cout ) {
			for ( int i = 0; i < (int)constraints.size(); ++i ) {
				out << i << ": " << constraints[i] << '\n';
			}
		}

		int GetBlockVariableValue( int x, int y ) {
			return cplex.getValue( blocksVars[base.index( x, y )] );
		}

		Base const & GetBase()       const {   return base;         }
		Base const & GetSolvedBase() const {   return solvedBase;   }

	private:
		struct RectVarsType {
			RectVarsType( IloBoolVarArray const & vars, int w, int h )
				: rectVars{ vars }, rectWidth{ w }, rectHeight{ h } {
			}
			IloBoolVarArray rectVars;
			int rectWidth, rectHeight;
		};

	private:
		IloEnv   env;
		IloModel model;
		IloCplex cplex;

		Base const & base;
		Base solvedBase;

		IloBoolVarArray        blocksVars;
		vector< RectVarsType > rectsVarsArray;

		vector< IloConstraint > constraints;

	private:
		void InitBaseVariables() {
			blocksVars = IloBoolVarArray{ env, base.size() };
			for ( int i = 0; i < base.size(); ++i ) {
				char buf[1024];
				snprintf( buf, 1024, "x_{%d,%d}", i%(base.width), i/base.width );
				blocksVars[i].setName( buf );

				BlockType blockType = base.blocks[i];
				if ( blockType == DoNotExistsBlock ) {
					AddConstraint( model, blocksVars[i] == 0 );
				}
			}
		}

		void InitRectVariables() {
			for ( int i = 1; i <= maxRectSize; ++i ) {
				for ( int j = 1; j <= maxRectSize; ++j ) {
					if ( i*j <= maxRectSize ) {
						rectsVarsArray.emplace_back( AddRectangle( env, model, blocksVars, i, j ), i, j );
					}
				}
			}
		}

		void InitOverlapConstraints() {
			for ( int i = 0; i < base.width; ++i ) {
				for ( int j = 0; j < base.height; ++j ) {
					IloExpr overlapExpr{ env };
					for ( int k = 0; k < (int)rectsVarsArray.size(); ++k ) {
						IloBoolVarArray const & rectVars = rectsVarsArray[k].rectVars;
						int w = rectsVarsArray[k].rectWidth;
						int h = rectsVarsArray[k].rectHeight;
						overlapExpr += ConstructInclusionConstraintsAndOverlapExpr( env, model, i, j
						                                                          , rectVars, w, h
						                                                          );
					}
					AddConstraint( model, overlapExpr <= 1 );
					AddConstraint( model, blocksVars[base.index( i, j )] <= overlapExpr );
				}
			}
		}

		void InitObjectiveFunction() {
			IloExpr objExpr{ env };
			for ( int i = 0; i < base.size(); ++i ) {
				objExpr += blocksVars[i];
			}
			model.add( IloObjective( env, objExpr, IloObjective::Maximize ) );
		}

		void SetCplexParameters( IloCplex & cplex ) {
			cplex.setParam( IloCplex::Param::Parallel, IloCplex::Opportunistic );
			cplex.setParam( IloCplex::Param::Threads, cplex.getNumCores() );
			cplex.setParam( IloCplex::Param::Emphasis::MIP, IloCplex::MIPEmphasisBestBound );
		}

		void FillSolvedBase() {
			solvedBase.width  = base.width;
			solvedBase.height = base.height;
			for ( int i = 0; i < base.size(); ++i ) {
				switch ( base.blocks[i] ) {
					case DoNotExistsBlock: {
						solvedBase.addBlock( DoNotExistsBlock );
					} break;
					case SelectedRegularBlock: // fall through
					case RegularBlock: {
						if ( cplex.getValue( blocksVars[i] ) == 0 ) {
							solvedBase.addBlock( RegularBlock );
						} else {
							solvedBase.addBlock( SelectedRegularBlock );
						}
					} break;
					case SelectedHamBlock: // fall through
					case HamBlock: {
						if ( cplex.getValue( blocksVars[i] ) == 0 ) {
							solvedBase.addBlock( HamBlock );
						} else {
							solvedBase.addBlock( SelectedHamBlock );
						}
					} break;
				}
			}
		}

		void AddConstraint( IloModel & model, IloConstraint const & constraint ) {
			model.add            ( constraint );
			constraints.push_back( constraint );
		}

		IloBoolVar AddMultConstraints( IloEnv & env, IloModel & model, IloBoolVar const & v0, IloBoolVar const & v1 ) {
			IloBoolVar res{ env };
			AddConstraint( model, res <= v0 );
			AddConstraint( model, res <= v1 );
			AddConstraint( model, res >= v0 + v1 - 1 );
			return res;
		}

		IloBoolVar AddProductConstraints( IloEnv & env, IloModel & model
		                                , IloBoolVarArray const & baseVars, int x, int y
		                                , int w, int h
		                                ) {
			queue< IloBoolVar > mulQueue;
			for ( int i = 0; i < w; ++i ) {
				for ( int j = 0; j < h; ++j ) {
					mulQueue.push( baseVars[x+i + base.width*(y+j)] );
				}
			}

			while ( mulQueue.size() > 1 ) {
				IloBoolVar v0 = mulQueue.front();
				mulQueue.pop();
				IloBoolVar v1 = mulQueue.front();
				mulQueue.pop();

				IloBoolVar mult = AddMultConstraints( env, model, v0, v1 );
				mulQueue.push( mult );
			}

			IloBoolVar last = mulQueue.front();
			mulQueue.pop();
			return last;
		}

		bool IsValidRect( int x, int y, int w, int h ) const {
			int sum = 0;
			for ( int i = 0; i < w; ++i ) {
				for ( int j = 0; j < h; ++j ) {
					BlockType blockType = base.access( x+i, y+j );
					if      ( blockType == DoNotExistsBlock ) {   return false;   }
					else if ( blockType == HamBlock         ) {   ++sum;          }
				}
			}
			return sum >= minNumCross;
		}

		IloBoolVarArray AddRectangle( IloEnv & env, IloModel & model, IloBoolVarArray const & baseVars, int w, int h ) {
			int matrixW = base.width-w+1;
			int matrixH = base.height-h+1;
			int matrixSize = matrixW*matrixH;
			IloBoolVarArray rect{ env, matrixSize };
			for ( int i = 0; i < matrixSize; ++i ) {
				char buf[1024];
				snprintf( buf, 1024, "%dx%d_{%d,%d}", w, h, i % matrixW, i / matrixW );
				rect[i].setName( buf );
			}

			for ( int i = 0; i < matrixW; ++i ) {
				for ( int j = 0; j < matrixH; ++j ) {
					if ( !IsValidRect( i, j, w, h ) ) {
						AddConstraint( model, rect[i + matrixW*j] == 0 );
					} else {
						IloBoolVar product = AddProductConstraints( env, model, baseVars, i, j, w, h );
						AddConstraint( model, rect[i + matrixW*j] <= product );
					}
				}
			}
			return rect;
		}

		IloExpr ConstructInclusionConstraintsAndOverlapExpr( IloEnv & env, IloModel & model, int x, int y
		                                                   , IloBoolVarArray const & rects, int w, int h
		                                                   ) {
			IloExpr expr{ env };
			for ( int i = 0; i < w; ++i ) {
				for ( int j = 0; j < h; ++j ) {
					if ( x-i >= 0 && x-i < base.width-w+1 && y-j >= 0 && y-j < base.height-h+1 ) {
						AddConstraint( model, blocksVars[base.index( x, y )] >= rects[x-i + (base.width-w+1)*(y-j)] );
						expr += rects[x-i + (base.width-w+1)*(y-j)];
					}
				}
			}
			return expr;
		}

		void OutputSolution( ostream & o ) const {
			cout << "Input: \n";
			PrintBase( base );
			cout << "\nSolution: \n";
			PrintBase( solvedBase );
			cout << "Score: " << GetBestScore() << endl;
		}

		void PrintBase( Base const & base ) const {
			cout << '+';
			for ( int i = 0; i < base.width; ++i ) {   cout << '=';   }
			cout << "+\n";
			for ( int j = 0; j < base.height; ++j ) {
				cout << '|';
				for ( int i = 0; i < base.width; ++i ) {   cout << blockStrings[base.access( i, j )];   }
				cout << "|\n";
			}
			cout << '+';
			for ( int i = 0; i < base.width; ++i ) {   cout << '=';   }
			cout << "+\n";
		}
};

vector< BlockType > CopySubMatrix( Base const & base, int x, int y, int w, int h ) {
	vector< BlockType > res;
	for ( int j = y; j < y+h; ++j ) {
		for ( int i = x; i < x+w; ++i ) {
			res.push_back( base.access( i, j ) );
		}
	}
/*
	for ( int i = 0; i < res.size(); ++i ) {
		if ( i % baseW == 0 ) cout << "\n";
		switch ( initialBase.first[i] ) {
			case REGULAR_TYPE: cout << "-"; break;
			case HAM_TYPE: cout << "o"; break;
			case DO_NOT_EXISTS_TYPE: cout << "."; break;
		}
	}
*/
	return res;
}

vector< Base > DivideBase( Base const & base, int numX, int numY ) {
	int w = base.width;
	int h = base.height;
	int wD = w / numX;
	int hD = h / numY;

	vector< Base > bases;
	for ( int i = 0; i < numX; ++i ) {
		for ( int j = 0; j < numY; ++j ) {
			int curX = i*wD, curY = j*hD;
			int curW, curH;
			if ( i == numX-1 ) {
				if ( j == numY-1 )     { curW = w-(i*wD);   curH = h-(j*hD);   }
				else                   { curW = w-(i*wD);   curH =      hD;    }
			} else if ( j == numY-1 ) { curW =      wD;    curH = h-(j*hD);   }
			else                      { curW =      wD;    curH =      hD;    }
			bases.emplace_back( CopySubMatrix( base, curX, curY, curW, curH ), curW, curH );
		}
	}
	return bases;
}

int main( int argc, char * argv[] ) {
	if ( argc != 5 ) {
		cerr << "Usage: " << argv[0] << " <file> <outfile> [<num subdivide x> <num subdivide y>]" << endl;
		exit( 0 );
	}

	Base initialBase = LoadBaseFromFile( argv[1] );

	int numSubdivideX = (argc > 3 ? atoi( argv[3] ) : 1);
	int numSubdivideY = (argc > 4 ? atoi( argv[4] ) : 1);
	auto bases = DivideBase( initialBase, numSubdivideX, numSubdivideY );
	vector< SubProblem > subProblems;
	for ( auto const & base: bases ) {
		subProblems.emplace_back( base );
	}
	int numIter = 0, sumScores = 0;
	for ( auto & problem: subProblems ) {
		cout << "===== Iteration number: " << ++numIter << " / " << bases.size() << " ======" << endl;
		problem.Solve();
		sumScores += problem.GetBestScore();
		cout << "===== Cumulative score: " << sumScores << " ======" << endl;
		cout << "===== Potential score: "  << subProblems.size()/(double)numIter*sumScores << " ======\n" << endl;
	}
	cout << "Total score w/ Do Not Exists: " << sumScores + initialBase.numDoNotExists << endl;

	return 0;
}
