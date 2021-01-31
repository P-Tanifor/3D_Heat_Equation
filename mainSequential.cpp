#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cmath>
using namespace std::chrono;

using namespace std;
vector<vector<vector<double> > > matrix3D,previousMatrix3D;
ifstream fin("ADI_input.txt");
int x, y, z;
int u0, k, N, M, Q;
double precision, dt;
double lamdaX , lamdaY , lamdaZ ;
double deltaX, deltaY, deltaZ;



void readData() {
    fin >> x >> y >> z ;
    fin >>u0 >> k;
    fin >> precision;
    fin >> N >>M >> Q;
    fin >> dt;
    deltaX = 1.0 * x / (N-1);
    deltaY = 1.0 * y / (M-1);
    deltaZ = 1.0 * z / (Q-1);
    lamdaX =1.0 / 3 * dt/(deltaX * deltaX);
    lamdaY =1.0 / 3 * dt/(deltaY * deltaY);
    lamdaZ =1.0 / 3 * dt/(deltaZ * deltaZ);


    cout << "x = " << x<< endl;
    cout << "y = " << y<< endl;
    cout << "z = " << z<< endl;
    cout << "u0 = " <<u0 <<endl;
    cout << "k = " <<k <<endl;
    cout << "precision = " << precision<<endl;
    cout << "N = " << N <<endl;
    cout << "M = " << M <<endl;
    cout << "Q = " << Q <<endl;
    cout << "dt = " << dt <<endl;
}

bool is_onborder(int i, int j, int l){

    return i == 0 || i==N-1 || j==0 || j==M-1 || l== 0 || l== Q-1;
}

void thomas(vector<double>& r, double lamda){
    int i, LEN;
    LEN =r.size();
    vector<double> a(LEN), b(LEN), c(LEN);
    b[0] = 1;
    b[LEN-1] = 1;
    for(i=1; i<= LEN-2; i++){
        a[i] = -lamda;
        b[i] = 1 + 2*lamda;
        c[i] = -lamda;
    }
    for(i=1; i<=LEN-2; i++){
        b[i] -= a[i] * c[i-1];
        r[i] -= a[i] * r[i-1];
        a[i] = 0;

        c[i] /= b[i];
        r[i] /= b[i];
        b[i] = 1;

    }
    /* back substitution */
    for(i=LEN-2; i>0; i--){
        r[i] -= r[i+1] * c[i];
        c[i] -= b[i+1] * c[i];
        //r[i] -= r[i+1] * c[i];
    }


}

void printMatrix(){
    int i, j, l;
    for (i=0; i<N; i++){
        for (j=0; j<M; j++){
            for (l=0; l<Q; l++){
                cout << matrix3D[i][j][l] << " ";
            }
            cout <<endl;
        }
        cout << endl;
    }
    cout <<"................."<<endl;

}
//void printVector(vector<double>& v){
//    for (auto x : v){
//        cout << " " << x ;
//    }
//    cout << endl;
//}

void copyMatrix(){
    int i, j, l;
    for (i=0; i<N; i++){
        for (j=0; j<M; j++){
            for (l=0; l<Q; l++){
                previousMatrix3D[i][j][l] = matrix3D[i][j][l];
            }
        }
    }

}

double error(){
    double result;
    double sumOfSquaredDiff = 0.0;
    int i, j, l;
    for (i=0; i<N; i++){
        for(j=0; j<M; j++){
            for(l=0; l<Q; l++){
                double diff;
                diff = matrix3D[i][j][l] - previousMatrix3D[i][j][l];
                double squaredDiff;
                squaredDiff = diff * diff;
                sumOfSquaredDiff += squaredDiff;
            }
        }
    }
    result = sqrt(sumOfSquaredDiff / (N * M * Q));
    return result;
}

int main() {

    auto start = high_resolution_clock::now();
    readData();
    int i,j,l;
    matrix3D.resize(N);
    previousMatrix3D.resize(N);
    for (i=0; i< N; i++){
        matrix3D[i].resize(M);
        previousMatrix3D[i].resize(M);
        for (j=0; j<M; j++){
            matrix3D[i][j].resize(Q);
            previousMatrix3D[i][j].resize(Q);
        }
    }
    for (i=0; i<N; i++){
        for (j=0; j<M; j++){
            for (l=0; l<Q; l++){
                matrix3D[i][j][l] = k;
                if (is_onborder(i,j,l)) {
                    matrix3D[i][j][l] = u0;
                }
            }
        }
    }
    //printMatrix();
    copyMatrix();
    int step = 1;
    while(true){
        cout << "step = " <<step <<endl;
        //fix directions Y and Z, then iterate over X
        for(j=1; j<M-1; j++){
            for(l=1; l<Q-1; l++){
                // we will copy the vector we will work on. matrix3D[0...N-1][j][l]
                vector<double> v;
                v.resize(N);
                for(i=0; i<N; i++){
                    v[i] =matrix3D[i][j][l];
                }
                //apply the thomas algorithm on it
                thomas(v, lamdaX);
                /* cpy results into new matrix */
                for(i=0; i<N; i++){
                    matrix3D[i][j][l] = v[i];
                }
            }
        }

        /* fix directions X and Z, then iterate over Y.*/
        for(i=1; i<N-1; i++){
            for(l=1; l<Q-1; l++){
                /* we will copy the vector we will work on. matrix3D[i][0...M-1][l]*/
                vector<double> v;
                v.resize(M);
                for(j=0; j<M; j++){
                    v[j] = matrix3D[i][j][l];
                }
                /* we apply the thomas method for it*/
                thomas(v, lamdaY);
                //cpy the results into new matrix
                for(j=0; j<M; j++){
                    matrix3D[i][j][l] = v[j];
                }
            }
        }
        /* we fix direction X and Y, then iterate over Z */
        for(i=1; i<N-1; i++){
            for(j=1; j<M-1; j++){
                //cpy the vector we will work on --> matrix3D[i][j][0...l-1]
                vector<double> v;
                v.resize(Q);
                for(l=0; l<Q; l++){
                    v[l] = matrix3D[i][j][l];
                }
                // apply the thomas method for it
                thomas(v,lamdaZ);
                // cpy the result into the new matrix
                for(l=0; l<Q; l++){
                    matrix3D[i][j][l] = v[l];
                }
            }
        }
        /** compare current matrix with the previous matrix, if the difference between them is
         * less than the precision, then stop the algorithm.
         **/
        double err = error();
        cout << "error at step "<<step <<" is " << err <<endl;
        if(err < precision){
            cout << "Finished with error "<< err << "!"<<endl;
            break;
        }
        copyMatrix();


        step ++;

    }
    printMatrix();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    cout << "The algorithm took " << 1.0 * duration.count() / 1000 << " seconds to execute!" << endl;



    return 0;
}


