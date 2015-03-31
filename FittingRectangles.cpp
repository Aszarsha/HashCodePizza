#define RUN_ME_WITH_SH /*
g++ -Wall -std=c++1y -O3 -DIL_STD -DVERBOSE                            \
	-o FittingRectangles FittingRectangles.cpp                          \
	-I/opt/ibm/ILOG/CPLEX_Studio1261/cplex/include/                      \
	-I/opt/ibm/ILOG/CPLEX_Studio1261/concert/include/                    \
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
#include <ilcplex/ilocplexi.h>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <queue>

using namespace std;

int minNumCross = 3;
int maxRectSize = 12;

vector< IloConstraint > constraints;

int numCross;
int numBlank;

typedef pair< vector< int >, pair< int, int > > Base;

Base initialBase;

int numNodesX, numNodesY, numNodes;
vector< int > const * baseSet;

Base ConstructBase( vector< int > const & base, int width, int height ) {
	return make_pair( base, make_pair( width, height ) );
}

void LoadBase( char const * fileName ) {
	FILE * fin = fopen( fileName, "rb" );
	int a, b, c, d;
	fscanf( fin, "%d %d %d %d\n", &a, &b, &c, &d );
	char rc;
	for ( int i = 0; i < (b+1)*a; ++i ) {
		fscanf( fin, "%c", &rc );
		if ( rc == '\n' ) continue;
		if ( rc == 'H' ) {   initialBase.first.push_back( 1 );   ++numCross;   }
		else             {   initialBase.first.push_back( 0 );   ++numBlank;   }
	}
	initialBase.second.first  = b;
	initialBase.second.second = a;
	minNumCross = c;
	maxRectSize = d;

	printf( "Blank: %d, total: %d, ratio: %lf\n", numCross, numCross+numBlank, numCross/(double)(numCross+numBlank) );
}

vector< int > CopySubMatrix( vector< int > const & base, int baseW, int x, int y, int w, int h ) {
	vector< int > res;
	for ( int i = x; i < x+w; ++i ) {
		for ( int j = y; j < y+h; ++j ) {
			res.push_back( base[i + baseW*j] );
		}
	}
	return res;
}

vector< Base > DivideBase( Base const & baseSet, int numX, int numY ) {
	int w = baseSet.second.first;
	int h = baseSet.second.second;
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
			bases.push_back( ConstructBase( CopySubMatrix( baseSet.first, w, curX, curY, curW, curH ), curW, curH ) );
		}
	}
	return bases;
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
			mulQueue.push( baseVars[x+i + numNodesX*(y+j)] );
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

int ComputeSum( int x, int y, int w, int h ) {
	int sum = 0;
	for ( int i = 0; i < w; ++i ) {
		for ( int j = 0; j < h; ++j ) {
			if ( (*baseSet)[x+i + numNodesX*(y+j)] ) {
				++sum;
			}
		}
	}
	return sum;
}

IloBoolVarArray AddRectangle( IloEnv & env, IloModel & model, IloBoolVarArray const & baseVars, int w, int h ) {
	cout << " " << w << "x" << h;
	int matrixW = (numNodesX-w+1);
	int matrixH = (numNodesY-h+1);
	int matrixSize = matrixW*matrixH;
	IloBoolVarArray rect{ env, matrixSize };
	for ( int i = 0; i < matrixSize; ++i ) {
		char buf[1024];
		snprintf( buf, 1024, "%dx%d_{%d,%d}", w, h, i % matrixW, i / matrixW );
		rect[i].setName( buf );
	}

	for ( int i = 0; i < matrixW; ++i ) {
		for ( int j = 0; j < matrixH; ++j ) {
			IloBoolVar product = AddProductConstraints( env, model, baseVars, i, j, w, h );
			IloConstraint c0 = rect[i + matrixW*j] <= product;
			model.add            ( c0 );
			constraints.push_back( c0 );

			IloConstraint c1 = rect[i + matrixW*j] <= (ComputeSum( i, j, w, h ) >= minNumCross ? 1 : 0);
			model.add            ( c1 );
			constraints.push_back( c1 );
		}
	}

	return rect;
}

IloExpr ConstructOverlapExpr( IloEnv & env, IloModel & model, int x, int y
                            , IloBoolVarArray const & v0, int w0, int h0
                            ) {
	IloExpr expr{ env };
	for ( int i = 0; i < w0; ++i ) {
		for ( int j = 0; j < h0; ++j ) {
			if ( x-i >= 0 && x-i < numNodesX-w0+1 && y-j >= 0 && y-j < numNodesY-h0+1 ) {
				expr += v0[x-i + (numNodesX-w0+1)*(y-j)];
			}
		}
	}
	return expr;
}

IloExpr ConstructInclusionConstraintsAndExpr( IloEnv & env, IloModel & model
                                            , IloBoolVarArray const & baseVars, int x, int y
                                            , IloBoolVarArray const & rects, int w, int h
                                            ) {
	IloExpr expr{ env };
	for ( int i = 0; i < w; ++i ) {
		for ( int j = 0; j < h; ++j ) {
			if ( x-i >= 0 && x-i < numNodesX-w+1 && y-j >= 0 && y-j < numNodesY-h+1 ) {
				model.add            ( baseVars[x + numNodesX*y] >= rects[x-i + (numNodesX-w+1)*(y-j)] );
				constraints.push_back( baseVars[x + numNodesX*y] >= rects[x-i + (numNodesX-w+1)*(y-j)] );
				expr += rects[x-i + (numNodesX-w+1)*(y-j)];
			}
		}
	}
	return expr;
}

typedef pair< IloBoolVarArray, pair< int, int > > RectVarsType;
typedef vector< RectVarsType > RectVarsArrayType;

RectVarsType ConstructVarsArrays( IloEnv & env, IloModel & model, IloBoolVarArray const & baseVars, int w, int h ) {
	auto res = make_pair( AddRectangle( env, model, baseVars, w, h ), make_pair( w, h ) );
	return res;
}


RectVarsArrayType AddAllRectangles( IloEnv & env, IloModel & model, IloBoolVarArray const & baseVars ) {
	RectVarsArrayType varsArrays;
	cout << "Adding rectangles";
#if 1
	for ( int i = 1; i <= maxRectSize; ++i ) {
		for ( int j = 1; j <= maxRectSize; ++j ) {
			if ( i*j <= maxRectSize ) {
				varsArrays.push_back( ConstructVarsArrays( env, model, baseVars, i, j ) );
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
	cout << endl;
	return varsArrays;
}

int sumScores;

void Solve() {
	IloEnv   env;   // Default init
	IloModel model{ env };
	IloCplex cplex{ model };

	IloBoolVarArray baseVars{ env, numNodes };
	for ( int i = 0; i < numNodes; ++i ) {
		char buf[1024];
		snprintf( buf, 1024, "x_{%d,%d}", i%numNodesX, i/numNodesX );
		baseVars[i].setName( buf );
	}

	RectVarsArrayType varsArray = AddAllRectangles( env, model, baseVars );

	for ( int i = 0; i < numNodesX; ++i ) {
		for ( int j = 0; j < numNodesY; ++j ) {
			IloExpr expr( env );
			for ( int k = 0; k < (int)varsArray.size(); ++k ) {
				IloBoolVarArray const & vars = varsArray[k].first;
				int w = varsArray[k].second.first;
				int h = varsArray[k].second.second;
				expr += ConstructOverlapExpr( env, model, i, j, vars, w, h );
			}
			model.add            ( expr <= 1 );
			constraints.push_back( expr <= 1 );
		}
	}

	for ( int i = 0; i < numNodesX; ++i ) {
		for ( int j = 0; j < numNodesY; ++j ) {
			IloExpr expr{ env };
			for ( int k = 0; k < (int)varsArray.size(); ++k ) {
				IloBoolVarArray const & vars = varsArray[k].first;
				int w = varsArray[k].second.first;
				int h = varsArray[k].second.second;
				expr += ConstructInclusionConstraintsAndExpr( env, model, baseVars, i, j, vars, w, h );
			}
			model.add            ( baseVars[i + numNodesX*j] <= expr );
			constraints.push_back( baseVars[i + numNodesX*j] <= expr );
		}
	}

	IloExpr objExpr{ env };
	for ( int i = 0; i < numNodes; ++i ) {
		objExpr += baseVars[i];
	}
	model.add( IloObjective( env, objExpr, IloObjective::Maximize ) );

#if 0
	for ( int i = 0; i < (int)constraints.size(); ++i ) {
		cout << i << ": " << constraints[i] << '\n';
	}
#endif

	if ( !cplex.solve() ) {
		if ( cplex.getStatus() != IloAlgorithm::Infeasible ) {
			clog << "Optimization error, CPLEX status code: " << cplex.getStatus() << endl;
		} else {
			clog << "Infeasible, CPLEX status code: " << cplex.getStatus() << endl;
		}
	} else {
		cout << "Input: \n";
		for ( int j = 0; j < numNodesY; ++j ) {
			for ( int i = 0; i < numNodesX; ++i ) {
				cout << ((*baseSet)[i + numNodesX*j] == 0 ? "-" : "•");
			}
			cout << endl;
		}
		cout << "Solution: \n";
		for ( int j = 0; j < numNodesY; ++j ) {
			for ( int i = 0; i < numNodesX; ++i ) {
				if ( cplex.getValue( baseVars[i + numNodesX*j] ) == 0 ) {
					cout << ((*baseSet)[i + numNodesX*j] == 1 ? "•" : "-");
				} else {
					cout << ((*baseSet)[i + numNodesX*j] == 1 ? "o" : "x");
				}
			}
			cout << endl;
		}

		cout << "Best objective: " << cplex.getBestObjValue() << endl;
		sumScores += cplex.getBestObjValue();
	}
}

int main( int argc, char * argv[] ) {
	if ( argc != 4 ) {
		cerr << "Usage: " << argv[0] << " <file> <num subdivide x> <num subdivide y>" << endl;
		exit( 0 );
	}

	LoadBase( argv[1] );

	auto bases = DivideBase( initialBase, atoi( argv[2] ), atoi( argv[3] ) );
	int numIter = 0;
	for ( auto const & base: bases ) {
		cout << "===== Iteration number: " << ++numIter << " / " << bases.size() << " ======" << endl;
		numNodesX = base.second.first;
		numNodesY = base.second.second;
		numNodes = numNodesX * numNodesY;
		baseSet = &base.first;
		Solve();
		cout << "===== Cumulative score: " << sumScores << " ======" << endl;
		cout << "===== Potential score: "  << bases.size()/(double)numIter * sumScores << " ======\n" << endl;
	}

	cout << "Total score: " << sumScores << endl;

	return 0;
}
