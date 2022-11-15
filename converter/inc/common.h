#ifndef COMM_H
#define COMM_H

/* dimention */ 
#define DIM 2

/* droplet number */
#define DROPLET_NUMBER 2

/* grid number */ 
#define Nl 18   // 液相球座標　r方向
#define Ngl 26  // 液滴大側　気相球座標　r方向 
#define Ngs 26  // 液滴小側　気相球座標　r方向
#define NQl 26  // 液滴大側　気相球座標　θ方向
#define NQs 26  // 液滴小側　気相球座標　θ方向
#define NZ 350  // 気相円筒座標　z方向
#define NX 150  // 気相円筒座標　r方向

/* donor point */   // 球座標系→円筒座標系への補間位置を決める．
#define DONOR_L 15
#define DONOR_S 15

/* initial droplet size */
#define Rsl0 50e-5  // [m]
#define Rss0 30e-5  // [m]

/* coordinate area */
#define rlmax 180e-5 // 液滴小側　気相球座標系　範囲
#define rsmax 180e-5 // 液滴大側　気相球座標系　範囲
#define Lmax 3500e-5 // 気相円筒座標系　z方向　範囲
#define Xmax 1500e-5 // 気相円筒座標系　r方向　範囲

/* CFL number */
#define MAXCRN 0.1 // 最大クーラン数
#define MAXDIF 0.40  // 最大拡散数 

/* output setting */
// #define ORDER 1000000.0
#define ORDER 1000.0
// #define OUTPUT_INTERVAL 1.0e-4
#define OUTPUT_INTERVAL 3.0e-3

/* iteration number */
#define Kmax 5000000
#define Mmax 30000
#define erefZ 1.0e-8
#define erefT 1.0e-5
#define erefY 1.0e-7
#define dtdef 5.0e-5
#define dtdefini 1.0e-7
#define dtinitstep 100
#define dtchagestep 2000
#define tmax 15

/* iteration number */
// #define Kmax 5000000
// #define Mmax 30000
// #define erefZ 1.0e-8
// #define erefT 1.0e-5
// #define erefY 1.0e-7
// #define dtdef 5.0e-6
// #define dtdefini 1.0e-8
// #define dtinitstep 100
// #define dtchagestep 2000
// #define tmax 10

/* point implicit method */
#define _Dif_implicit_flg 0

typedef int Droplet;

enum {
    SMALL,
    LARGE,
};

#endif