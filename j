13a14
> 
64c65
< void GeneralizedEigenvalueProblem(long int N, double * A, double * B, double * c, double * eig, bool * energy_is_real){
---
> void GeneralizedEigenvalueProblem(long int N, double * A, double * B, double * c, double * eig){
81,99c82,102
< //    outfile->Printf("\n");
< //    outfile->Printf("Before solving eq, the matrix A is\n");
< //    for(int i = 0; i<N;i++){
< //       for(int j = 0; j<N;j++){
< //       Al->pointer()[i][j] = A[i*N+j];
< //       }
< //    }
< //    Al->print();
< //     outfile->Printf("\n");
< //    std::shared_ptr<Matrix> Bl (new Matrix(N,N));
< //    outfile->Printf("\n");
< //    outfile->Printf("Before solving eq, the matrix B is\n");
< //    for(int i = 0; i<N;i++){
< //       for(int j = 0; j<N;j++){
< //       Bl->pointer()[i][j] = B[i*N+j];
< //       }
< //    }
< //    Bl->print();
< //     outfile->Printf("\n");
---
>     outfile->Printf("\n");
>     outfile->Printf("Before solving eq, the matrix A is\n");
>     for(int i = 0; i<N;i++){
>        for(int j = 0; j<N;j++){
>        Al->pointer()[i][j] = A[i*N+j];
>        }
>     }
>     Al->print();
>      outfile->Printf("\n");
> 
>     std::shared_ptr<Matrix> Bl (new Matrix(N,N));
>     outfile->Printf("\n");
>     outfile->Printf("Before solving eq, the matrix B is\n");
>     for(int i = 0; i<N;i++){
>        for(int j = 0; j<N;j++){
>        Bl->pointer()[i][j] = B[i*N+j];
>        }
>     }
>     Bl->print();
>      outfile->Printf("\n");
> 
115,117c118
<        // eig[i] = ALPHAR[i] / BETA[i];
<        // alphai->pointer()[i]=ALPHAI[i];
<       //  beta->pointer()[i]=BETA[i];
---
>         eig[i] = ALPHAR[i] / BETA[i];
119,127d119
<        if ( fabs(ALPHAI[i]/BETA[i]) > 1e-6 ) {
<            //eig[i] = 0.0;
<            eig[i] = ALPHAR[i] / BETA[i];
<            outfile->Printf("<<<warning>>> excitation energy is complex: %20.12lf + %20.12lf I\n",eig[i], ALPHAI[i]/BETA[i]);
<            energy_is_real[i] = false;
<        }else {
<            eig[i] = ALPHAR[i] / BETA[i];
<            energy_is_real[i] = true;
<        }
138,156c130,148
< //    outfile->Printf("\n");
< //    outfile->Printf("After solving eigenvalue problem, the matrix A is\n");
< //    std::shared_ptr<Matrix> An (new Matrix(N,N));
< //    for(int i = 0; i<N;i++){
< //       for(int j = 0; j<N;j++){
< //       An->pointer()[i][j] = A[i*N+j];
< //       }
< //    }
< //    An->print();
< //
< //    outfile->Printf("\n");
< //    outfile->Printf("After solving eigenvalue problem, the matrix B is\n");
< //    std::shared_ptr<Matrix> Bn (new Matrix(N,N));
< //    for(int i = 0; i<N;i++){
< //       for(int j = 0; j<N;j++){
< //       Bn->pointer()[i][j] = B[i*N+j];
< //       }
< //    }
< //    Bn->print();
---
>     outfile->Printf("\n");
>     outfile->Printf("After solving eigenvalue problem, the matrix A is\n");
>     std::shared_ptr<Matrix> An (new Matrix(N,N));
>     for(int i = 0; i<N;i++){
>        for(int j = 0; j<N;j++){
>        An->pointer()[i][j] = A[i*N+j];
>        }
>     }
>     An->print();
> 
>     outfile->Printf("\n");
>     outfile->Printf("After solving eigenvalue problem, the matrix B is\n");
>     std::shared_ptr<Matrix> Bn (new Matrix(N,N));
>     for(int i = 0; i<N;i++){
>        for(int j = 0; j<N;j++){
>        Bn->pointer()[i][j] = B[i*N+j];
>        }
>     }
>     Bn->print();
187c179,190
< 
---
>     std::shared_ptr<Matrix> Bl (new Matrix(N,N));
> //    outfile->Printf("\n");
> //    outfile->Printf("Before solving eq, the matrix B is\n");
> //    for(int i = 0; i<N;i++){
> //       for(int j = 0; j<N;j++){
> //       //	  double dum = 0.0;
> //         // dum += B[i*N+j];
> //       Bl->pointer()[i][j] = B[i*N+j];
> //       }
> //    }
> //    Bl->print();
> //     outfile->Printf("\n");
189a193,208
> //    outfile->Printf("Eigenvalues\n");
> //    for(int i = 0; i<N;i++){
> //        outfile->Printf("%20.12lf\n",eig[i]);
> //    }
> //     outfile->Printf("\n");
> //        outfile->Printf("After solving eigenvalue problem, the matrix B is\n");
> // std::shared_ptr<Matrix> Bn (new Matrix(N,N));
> //    for(int i = 0; i<N;i++){
> //       for(int j = 0; j<N;j++){
> //         // double dum = 0.0;
> //         // dum += B[i*N+j];
> //       Bn->pointer()[i][j] = B[i*N+j];
> //       }
> //    }
> //    Bn->print();
>  
