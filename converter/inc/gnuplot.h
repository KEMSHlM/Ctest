#ifndef PLOT_H
#define PLOT_H

#include "common.h"

/*　描画用静止時間 */
#define SLEEP_TIME 0.05
/* データ出力間隔　
*   1で全てを出力(3次元では見にくい)
*/
#define OUP_INTERVAL 1

typedef int Coordinate;
typedef int Viewmode;

enum {
    SCALAR_SPHERICAL_S, // 液滴(小)球座標系のスカラー
    SCALAR_SPHERICAL_L, //　液滴(大)球座標系のスカラー
    SCALAR_CYLINDRICAL, //  円筒座標系のスカラー
    VECTOR_SPHERICAL_QS, //　液滴(小)球座標系の周方向ベクトル
    VECTOR_SPHERICAL_RS, //　液滴(小)球座標系の半径方向ベクトル
    VECTOR_SPHERICAL_QL, // 液滴(大)球座標系の周方向ベクトル　
    VECTOR_SPHERICAL_RL, // 液滴(大)球座標系の半径方向ベクトル
    VECTOR_CYLINDRICAL_Z, //　円筒座標系のz方向ベクトル
    VECTOR_CYLINDRICAL_X, //　円筒座標系の半径方向ベクトル
    VIEW_2D_COLOR = 0, //　色付き２次元マップでの可視化
    VIEW_2D_COLOR_LOG, //　色付きログスケール２次元マップでの可視化
    VIEW_2D_GREY, //　白黒２次元マップでの可視化
    VIEW_3D, //　3次元可視化
    VIEW_3D_COLOR, //　3次元可視化カラーマップ付き
    VIEW_3D_COLOR_LOG, // ログスケール3次元可視化カラーマップ付き
    VIEW_3D_GREY, // 3次元可視化白黒マップ付き
};

/* 初期化 */
void plot_init();
/*  可視化
*   p           ->  可視化したいパラメーター   
*   coord       ->  パラメーターpの種類（スカラーorベクトル，球座標系or円筒座標系，etc．）
*   view        ->  可視化モード選択
*   file_name   ->  ファイル名
*   k           ->  ファイル名識別番号
*/
void plot(double **p, Coordinate coord, Viewmode view, char *file_name, int k);

#endif