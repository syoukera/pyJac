#include <stdio.h>
#include <stdlib.h>
#include "jacob.h"

#define N_SPECIES 33
#define N_STATES 34

// LAPACK関数のプロトタイプ宣言
extern void dgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda, 
    double* wr, double* wi, double* vl, int* ldvl, 
    double* vr, int* ldvr, double* work, int* lwork, int* info);

// static const int n_species = 10;
// static const int n_reactions = 40;

int main() {
    const double t = 0.0;
    const double pres = 101325.0;
    
    // mass fraction H2:O2 phi 1.0
    double y[N_STATES] = {
       1.00000000e+03, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       1.96270854e-01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 6.54236179e-02, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 7.38305530e-01
    };
        
    double jac[N_SPECIES*N_SPECIES] = {0.0};
    double jac_lapack[N_STATES*N_STATES] = {0.0};
    
    // 固有値計算用の変数
    double wr[N_STATES], wi[N_STATES];  // 実部と虚部
    double vl[N_STATES * N_STATES];  // 左固有ベクトル（不要ならNULLでも可）
    double vr[N_STATES * N_STATES];  // 右固有ベクトル（不要ならNULLでも可）
    int lwork = 4 * N_STATES;  // 作業配列のサイズ
    double work[4 * N_STATES];  // 作業配列
    int info;
    
    char jobvl = 'N';  // 左固有ベクトルは計算しない
    char jobvr = 'N';  // 右固有ベクトルは計算しない
    int n = N_STATES;
    int lda = N_STATES;
    int ldvl = N_STATES;
    int ldvr = N_STATES;

    printf("Calling eval_jacob...\n");
    
    eval_jacob(t, pres, y, jac);

    printf("eval_jacob finished.\n");

    // 結果を表示
    printf("Output values: ");
    for (int i = 0; i < N_STATES; i++) {
        printf("%e ", jac[i + i*N_STATES]);
        // for (int j = 0; j < N_SPECIES; j++) {
        //     printf("%e ", jac[i + j*N_SPECIES]);
        // }
        // printf("\n");
    }
    printf("\n");
    
    // assign value to jac_lapack
    for (int i = 0; i < N_SPECIES; i++) {
        for (int j = 0; j < N_SPECIES; j++) {
            jac_lapack[i*N_STATES + j] = jac[i*N_SPECIES + j];
        }
    }

    // LAPACK の dgeev を呼び出し（ヤコビアンの固有値を計算）
    dgeev_(&jobvl, &jobvr, &n, jac_lapack, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info);
    
    if (info == 0) {
        printf("固有値（実部 虚部）:\n");
        for (int i = 0; i < N_STATES; i++) {
            printf("%e + %ei\n", wr[i], wi[i]);
        }
    } else {
        printf("固有値計算に失敗しました (info = %d)\n", info);
    }

    return 0;
}
