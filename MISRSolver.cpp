#define RUN_ME_WITH_SH /*
g++ -Wall -std=c++14 -g -DDEBUG -DIL_STD -DVERBOSE -o ${0%.*} ${0} \
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

// For now suppose those two are constants once the input file is loaded
int minNumCross = 3, maxRectSize = 12;

struct Block {
	public:
		enum Type { RegularType = 0
		          , HamType
		          , SelectedRegularType
		          , SelectedHamType
		          , DoNotExistsType
		          };

		static Block const Regular;
		static Block const SelectedRegular;
		static Block const Ham;
		static Block const SelectedHam;
		static Block const DoNotExists;

		Type const type;

		Block( Type t ): type{ t } {   }

		char const * GetString() {   return blockStrings[type];   }

	private:
		static char const * const blockStrings[];
};
char const * const Block::blockStrings[] = { "-", "+", "x", "o", "#" };
Block const Block::Regular         = Block{ Block::RegularType };
Block const Block::SelectedRegular = Block{ Block::SelectedRegularType };
Block const Block::Ham             = Block{ Block::HamType };
Block const Block::SelectedHam     = Block{ Block::SelectedHamType };
Block const Block::DoNotExists     = Block{ Block::DoNotExistsType };

struct Base {
	Base() {   }

	Base( vector< Block > && b, int w, int h )
		: blocks{ move( b ) }, width{ w }, height{ h } {
	}

	Base( vector< Base > const & bases, int w, int h )
		: width{ w }, height{ h } {
		int numBaseInWidth = 0, accWidth = 0;
		for ( int i = 0; i < width; ++i ) {
			if ( i >= accWidth + bases[numBaseInWidth].width ) {
				accWidth += bases[numBaseInWidth].width;
				++numBaseInWidth;
			}
		}
		++numBaseInWidth; // Count last one

		int accHeight = 0;
		int startingBase = 0;
		for ( int j = 0; j < height; ++j ) {
			if ( j >= accHeight + bases[startingBase].height ) {
				accHeight += bases[startingBase].height;
				startingBase += numBaseInWidth;
			}
			int y = j - accHeight;
			int accWidth = 0;
			int baseOffset = 0;
			for ( int i = 0; i < width; ++i ) {
				if ( i >= accWidth + bases[startingBase + baseOffset].width ) {
					accWidth += bases[startingBase + baseOffset].width;
					++baseOffset;
				}
				int x = i - accWidth;
				AddBlock( bases[startingBase + baseOffset].access( x, y ) );
			}
		}
	}

	void AddBlock( Block block ) {
		switch ( block.type ) {
			case Block::DoNotExistsType:     {   ++numDoNotExists;   break;   }
			case Block::SelectedRegularType: {   ++numRegular;       break;   }
			case Block::RegularType:         {   ++numRegular;       break;   }
			case Block::SelectedHamType:     {   ++numHam;           break;   }
			case Block::HamType:             {   ++numHam;           break;   }
		}
		blocks.push_back( block );
	}

	int size() const {
		return (int)blocks.size();
	}

	int index( int i, int j ) const {
		return i + width*j;
	}

	Block access( int i, int j ) const {
		return blocks[index( i, j )];
	}

	friend ostream & operator<<( ostream & o, Base const & b ) {
		o << '+';
		for ( int i = 0; i < b.width; ++i ) {   o << '=';   }
		o << "+\n";
		for ( int j = 0; j < b.height; ++j ) {
			o << '|';
			for ( int i = 0; i < b.width; ++i ) {   o << b.access( i, j ).GetString();   }
			o << "|\n";
		}
		o << '+';
		for ( int i = 0; i < b.width; ++i ) {   o << '=';   }
		o << "+\n";
		return o;
	}

	vector< Block > blocks;
	int width, height;
	int numRegular = 0;
	int numHam     = 0;
	int numDoNotExists = 0;
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
			case '\n': {                                           continue;    }
			case  '-': // fall through
			case  'T': {   base.AddBlock( Block::Regular         );   break;    }
			case  '+': // fall through
			case  'H': {   base.AddBlock( Block::Ham             );   break;    }
			case  'x': {   base.AddBlock( Block::SelectedRegular );   break;    }
			case  'o': {   base.AddBlock( Block::SelectedHam     );   break;    }
			case  '#': {   base.AddBlock( Block::DoNotExists     );   break;    }
			default:   {
				cerr << "Unknown character: '" << rc << "' in input file" << endl;
				exit( 1 );
			}
		}
	}
	fclose( fin );
	return base;
}

// Maximum Independent Set of Rectangles
class MISRSolver {
	public:
		MISRSolver( Base const & b )
			: env{}
			, model{ env }
			, cplex{ model }
			, base{ b } {

			InitBaseVariables();
			InitRectVariables();
			InitMIPStart();
			InitOverlapConstraints();
			InitObjectiveFunction();

			SetCplexParameters( cplex );
		}

		bool Solve() {
			if ( !cplex.solve() ) {
				if ( cplex.getStatus() != IloAlgorithm::Infeasible ) {
					clog << "Optimization error, CPLEX status code: " << cplex.getStatus() << endl;
				} else {
					clog << "Infeasible, CPLEX status code: " << cplex.getStatus() << endl;
				}
				return false;
			}
			FillSolvedBase();
			return true;
		}

		double GetBestScore() const {
			return cplex.getBestObjValue();
		}

		void PrintConstraints( ostream & o = cout ) {
			IloModel::Iterator it{ model };
			while ( it.ok() ) {
				if ( (*it).asConstraint().getImpl() ) {   o << (*it).asConstraint() << endl;   }
				++it;
			}
		}

		void PrintSolution( ostream & o ) const {
			o << "Input: \n" << base << "\nSolution: \n" << solvedBase << "Score: " << GetBestScore() << endl;
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

	private:
		void InitBaseVariables() {
			blocksVars = IloBoolVarArray{ env, base.size() };
			for ( int i = 0; i < base.size(); ++i ) {
				char buf[1024];
				snprintf( buf, 1024, "x_{%d,%d}", i%(base.width), i/base.width );
				blocksVars[i].setName( buf );

				if ( base.blocks[i].type == Block::DoNotExistsType ) {
					model.add( blocksVars[i] == 0 );
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

		void InitMIPStart() {
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
					model.add( overlapExpr <= 1 );
					model.add( blocksVars[base.index( i, j )] <= overlapExpr );
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
				switch ( base.blocks[i].type ) {
					case Block::DoNotExistsType: {
						solvedBase.AddBlock( Block::DoNotExists );
					} break;
					case Block::SelectedRegularType: // fall through
					case Block::RegularType: {
						if ( cplex.getValue( blocksVars[i] ) == 0 ) {
							solvedBase.AddBlock( Block::Regular );
						} else {
							solvedBase.AddBlock( Block::SelectedRegular );
						}
					} break;
					case Block::SelectedHamType: // fall through
					case Block::HamType: {
						if ( cplex.getValue( blocksVars[i] ) == 0 ) {
							solvedBase.AddBlock( Block::Ham );
						} else {
							solvedBase.AddBlock( Block::SelectedHam );
						}
					} break;
				}
			}
		}

		IloBoolVar AddMultConstraints( IloEnv & env, IloModel & model, IloBoolVar const & v0, IloBoolVar const & v1 ) {
			IloBoolVar res{ env };
			model.add( res <= v0 );
			model.add( res <= v1 );
			model.add( res >= v0 + v1 - 1 );
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
					Block::Type blockType = base.access( x+i, y+j ).type;
					if ( blockType == Block::DoNotExistsType ) {
						return false;
					} else if ( blockType == Block::HamType || blockType == Block::SelectedHamType ) {
						++sum;
					}
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
						model.add( rect[i + matrixW*j] == 0 );
					} else {
						IloBoolVar product = AddProductConstraints( env, model, baseVars, i, j, w, h );
						model.add( rect[i + matrixW*j] <= product );
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
						model.add( blocksVars[base.index( x, y )] >= rects[x-i + (base.width-w+1)*(y-j)] );
						expr += rects[x-i + (base.width-w+1)*(y-j)];
					}
				}
			}
			return expr;
		}
};

vector< Block > CopySubMatrix( Base const & base, int x, int y, int w, int h ) {
	vector< Block > res;
	for ( int j = y; j < y+h; ++j ) {
		for ( int i = x; i < x+w; ++i ) {
			res.push_back( base.access( i, j ) );
		}
	}
	return res;
}

vector< Base > DivideBase( Base const & base, int numX, int numY ) {
	int w = base.width;
	int h = base.height;
	int wD = w / numX;
	int hD = h / numY;

	vector< Base > bases;
	for ( int j = 0; j < numY; ++j ) {
		for ( int i = 0; i < numX; ++i ) {
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
	vector< MISRSolver > solvers;
	int numIter = 0, sumScores = 0;
	for ( auto const & base: bases ) {
		solvers.emplace_back( base );
		cout << "===== Iteration number: " << ++numIter << " / " << bases.size() << " ======" << endl;
		solvers.back().Solve();
		solvers.back().PrintSolution( cout );
		sumScores += solvers.back().GetBestScore();
		cout << "===== Cumulative score: " << sumScores << " ======" << endl;
		cout << "===== Potential score: "  << bases.size()/(double)numIter*sumScores << " ======\n" << endl;
	}
	cout << "Total score w/ Do Not Exists: " << sumScores + initialBase.numDoNotExists << endl;

	vector< Base > solvedBases;
	for ( auto const & solver: solvers ) {
		solvedBases.push_back( solver.GetSolvedBase() );
	}
	Base endBase{ solvedBases, initialBase.width, initialBase.height };
	cout << "Initial base: \n" << initialBase << "Solved base: \n" << endBase << "Score: " << sumScores << endl;

	return 0;
}
