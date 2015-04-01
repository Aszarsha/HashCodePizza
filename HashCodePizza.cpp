#define RUN_ME_WITH_SH /*
g++ -Wall -std=c++14 -Ofast -DIL_STD -DVERBOSE -o ${0%.*} ${0}          \
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

enum BlockType { REGULAR_TYPE, HAM_TYPE, DO_NOT_EXISTS_TYPE };

struct Base {
	Base() {   }
	Base( vector< BlockType > && b, int w, int h )
		: blocks{ move( b ) }, width{ w }, height{ h } {
	}

	vector< BlockType > blocks;
	int width, height;
	int numCross = 0;
	int numBlank = 0;
	int numDoNotExists = 0;

	int size() const {   return (int)blocks.size();   }
	int accessIndex( int i, int j ) const  {   return i + width*j;   }
	BlockType access( int i, int j ) const {   return blocks[accessIndex(i, j)];   }
};

Base LoadBaseFromFile( char const * fileName ) {
	FILE * fin = fopen( fileName, "rb" );

	Base base;
	fscanf( fin, "%d %d %d %d\n", &base.height, &base.width, &minNumCross, &maxRectSize );
	char rc;
	for ( int i = 0; i < (base.width+1)*base.height; ++i ) {
		fscanf( fin, "%c", &rc );
		if      ( rc == '\n' ) continue;
		else if ( rc == 'T'  ) {   base.blocks.push_back( REGULAR_TYPE       );   ++base.numBlank;         }
		if      ( rc == 'H'  ) {   base.blocks.push_back( HAM_TYPE           );   ++base.numCross;         }
		else if ( rc == '.'  ) {   base.blocks.push_back( DO_NOT_EXISTS_TYPE );   ++base.numDoNotExists;   }
	}
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
				if ( blockType == DO_NOT_EXISTS_TYPE ) {
					AddConstraint( model, blocksVars[i] == 0 );
				}
			}
		}

		void InitRectVariables() {
		#if 1
			for ( int i = 1; i <= maxRectSize; ++i ) {
				for ( int j = 1; j <= maxRectSize; ++j ) {
					if ( i*j <= maxRectSize ) {
						rectsVarsArray.emplace_back(
							AddRectangle( env, model, blocksVars, i, j ),
							i, j
						);
					}
				}
			}
		#else
			varsArrays.push_back( ConstructVarsArrays( env, model, baseVars, 1 , 12 ) );
			varsArrays.push_back( ConstructVarsArrays( env, model, baseVars, 12, 1  ) );
			varsArrays.push_back( ConstructVarsArrays( env, model, baseVars, 2 , 6  ) );
			varsArrays.push_back( ConstructVarsArrays( env, model, baseVars, 6 , 2  ) );
			varsArrays.push_back( ConstructVarsArrays( env, model, baseVars, 3 , 4  ) );
			varsArrays.push_back( ConstructVarsArrays( env, model, baseVars, 4 , 3  ) );
		#endif
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
					AddConstraint( model, blocksVars[base.accessIndex( i, j )] <= overlapExpr );
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
					if      ( blockType == DO_NOT_EXISTS_TYPE ) {   return false;   }
					else if ( blockType == HAM_TYPE           ) {   ++sum;          }
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
						AddConstraint( model, blocksVars[base.accessIndex( x, y )] >= rects[x-i + (base.width-w+1)*(y-j)] );
						expr += rects[x-i + (base.width-w+1)*(y-j)];
					}
				}
			}
			return expr;
		}

		void OutputSolution( ostream & o ) const {
			auto UnselectedGetString = [&]( int i, int j ) {
				BlockType blockType = base.access( i, j );
				switch ( blockType ) {
					case REGULAR_TYPE:       {   return "-";   }
					case HAM_TYPE:           {   return "•";   }
					case DO_NOT_EXISTS_TYPE: {   return " ";   }
				}
				return "ERROR";
			};

			auto SelectedGetString = [&]( int i, int j ) {
				BlockType blockType = base.access( i, j );
				switch ( blockType ) {
					case REGULAR_TYPE:       {   return "x";   }
					case HAM_TYPE:           {   return "o";   }
					case DO_NOT_EXISTS_TYPE: {   return " ";   }
				}
				return "ERROR";
			};

			cout << "Input: \n";
			PrintBase( [&]( int i, int j ) {
				return UnselectedGetString( i, j );
			});
			cout << "\nSolution: \n";
			PrintBase( [&]( int i, int j ) {
				if ( cplex.getValue( blocksVars[base.accessIndex( i, j )] ) == 0 ) {
					return UnselectedGetString( i, j );
				} else {
					return SelectedGetString( i, j );
				}
			});

			cout << "Best objective: " << GetBestScore() << endl;
		}

		template< typename GetElemenStringFunc >
		void PrintBase( GetElemenStringFunc && getElemenString ) const {
			cout << '+';
			for ( int i = 0; i < base.width; ++i ) {
				cout << '=';
			}
			cout << "+\n";
		
			for ( int j = 0; j < base.height; ++j ) {
				cout << '|';
				for ( int i = 0; i < base.width; ++i ) {
					cout << getElemenString( i, j );
				}
				cout << "|\n";
			}
		
			cout << '+';
			for ( int i = 0; i < base.width; ++i ) {
				cout << '=';
			}
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

	cout << "Total score: " << sumScores << endl;
	cout << "Total score w/ Do Not Exists: " << sumScores + initialBase.numDoNotExists << endl;

	return 0;
}
