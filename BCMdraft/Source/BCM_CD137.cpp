#include "Include/BCM_CD137.h"
#include <tuple>

const std::vector<std::string> interNames=
{"A_0","A_a","A_s","A_Ab","A_s_Ab","N_0","N_a","N_s","N_Ab","N_s_Ab","T_ns","T_0","T_s","T_bl","a","g","Ag"};

enum{A_0,A_a,A_s,A_Ab,A_s_Ab,N_0,N_a,N_s,N_Ab,N_s_Ab,T_ns,T_0,T_s,T_bl,a,g,Ag};

const std::vector<std::string> paramNames{"rA","rN","rT","rTsp",
                                          "gA0","gAa","gIA","aA0","aAa","aIA",
                                          "muA0","muAa","muIA","muaA","KamuA","kappaA",
                                          "rgA0","rgAa","raA0","rRLA0","k_AxA","k_AxN","k_AxT","k_AxAg",
                                          "k_AxAg_ag","k_AxAb","Kg_A0Aa","Ka_A0Aa","gN0","gNa","gIN","aN0",
                                          "aNa","aIN","muN0","muNa","muaN","KamuN","kappaN0","kappaNa","rgN0",
                                          "raN0","raNa","rRLN0","rRLNa","k_NxN","k_NxAb","k_N0","KA_N0_Na","Kg_N0_Na",
                                          "Ka_N0_Na","gT0","gTs","gIT","aT0","aTs","aIT","muT0","muTs","muIT","muaT",
                                          "KamuT","kappaT0","kappaTs","kappaIT","rgT0","rgTs","raT0","raTs","rRT0",
                                          "K_TxAb","tAp","mua","mug","muAg","kappa_Tym","maxCells","k_axAb","k_gxAb"};


enum {rA, rN, rT, rTsp, gA0, gAa, gIA, aA0, aAa, aIA,
      muA0, muAa, muIA, muaA,
      KamuA, kappaA,
      rgA0, rgAa, raA0, rRLA0, k_AxA, k_AxN, k_AxT, k_AxAg,
      k_AxAg_ag, k_AxAb, Kg_A0Aa, Ka_A0Aa, gN0, gNa, gIN, aN0,
      aNa, aIN, muN0, muNa, muaN, KamuN, kappaN0, kappaNa, rgN0,
      raN0, raNa, rRLN0, rRLNa, k_NxN, k_NxAb, k_N0NA_A, KA_N0_Na, Kg_N0_Na,
      Ka_N0_Na, gT0, gTs, gIT, aT0, aTs, aIT, muT0, muTs, muIT, muaT,
      KamuT, kappaT0, kappaTs, kappaIT, rgT0, rgTs, raT0, raTs, rRT0,
      K_TxAb, tAp, mua, mug, muAg, sTym, maxCells, k_axAb,k_gxAb};

const std::vector<std::string> observNames
{"a_obs","g_obs","rRLA_obs","rRLN_obs","rRT_obs","rg_A_obs","ra_A_obs","rg_N_obs","ra_N_obs",
  "rg_T_obs","ra_T_obs","prolifT_obs","apoptosisT_obs","nCells_obs"};

enum {a_obs,g_obs,rRLA_obs,rRLN_obs,rRT_obs,rg_A_obs,ra_A_obs,rg_N_obs,ra_N_obs,rg_T_obs,ra_T_obs,
      apopRateT_obs,prolifRate_obs,nCells_obs};

const std::vector<std::string> externNames{"nC0","Ag0","Ab_R","Ab_a","Ab_g","a_rec","g_rec"};

enum {nC0,Ag0,Ab_R,Ab_a,Ab_g,a_rec,g_rec};


std::vector<VariableInformation>  BCM_CD137::internalStateVariables()const
{
  return {{A_0,"A_0","nCells",LOG10},
{A_a,"A_a","nCells",LOG10},
{A_s,"A_s","nCells",LOG10},
{A_Ab,"A_Ab","nCells",LOG10},
{A_s_Ab,"A_s_Ab","nCells",LOG10},
{N_0,"N_0","nCells",LOG10},
{N_a,"N_a","nCells",LOG10},
{N_s,"N_s","nCells",LOG10},
{N_Ab,"N_Ab","nCells",LOG10},
{N_s_Ab,"N_s_Ab","nCells",LOG10},
{T_ns,"T_ns","nCells",LOG10},
{T_0,"T_0","nCells",LOG10},
{T_s,"T_s","nCells",LOG10},
{T_bl,"T_bl","nCells",LOG10},
{a,"a","ug/ml",LOG10},
{g,"g","ug/ml",LOG10},
{Ag,"Ag","ug/ml",LOG10}
};
}

std::vector<VariableInformation>  BCM_CD137::parameterVariables()const
{
  return {
  {rA,"rA","ncells/ncells",LOG10RATIO},
{ rN,"rN","ncells/ncells",LOG10RATIO},
{ rT,"rT","ncells/ncells",LOG10RATIO},
{ rTsp,"rTsp","ncells/ncells",LOG10RATIO},
{ gA0,"gA0","ug/ncells",LOG10},
{ gAa,"gAa","ug/ncells",LOG10},
{ gIA,"gIA","ug/ncells",LOG10},
{ aA0,"aA0","ug/ncells",LOG10},
{ aAa,"aAa","ug/ncells",LOG10},
{ aIA,"aIA","ug/ncells",LOG10},
{muA0,"muA0","1/h",LOG10},
{ muAa,"muAa","1/h",LOG10},
{ muIA,"muIA","1/h",LOG10},
{ muaA,"muaA","1/h",LOG10},
{KamuA,"KamuA","1/h",LOG10},
{ kappaA,"kappaA","1/h",LOG10},
{rgA0,"rgA0","1/h",LOG10RATIO},
 { rgAa,"rgAa","1/h",LOG10RATIO},
 { raA0,"raA0","1/h",LOG10RATIO},
 { rRLA0,"rRLA0","1/h",LOG10RATIO},
 { k_AxA,"k_AxA","1/h",LOG10},
 { k_AxN,"k_AxN","1/h",LOG10},
 { k_AxT,"k_AxT","1/h",LOG10},
 { k_AxAg,"k_AxAg","1/h",LOG10},
 {k_AxAg_ag,"k_AxAg_ag","1/h",LOG10},
 { k_AxAb,"k_AxAb","1/h",LOG10},
 { Kg_A0Aa,"Kg_A0Aa","1/h",LOG10},
 { Ka_A0Aa,"Ka_A0Aa","1/h",LOG10},
 { gN0,"gN0","1/h",LOG10},
 { gNa,"gNa","1/h",LOG10},
 { gIN,"gIN","1/h",LOG10},
 { aN0,"aN0","1/h",LOG10},
 {aNa,"aNa","1/h",LOG10},
 { aIN,"aIN","1/h",LOG10},
 { muN0,"muN0","1/h",LOG10},
 { muNa,"muNa","1/h",LOG10},
 { muaN,"muaN","1/h",LOG10},
 { KamuN,"KamuN","1/h",LOG10},
 { kappaN0,"kappaN0","1/h",LOG10},
 { kappaNa,"kappaNa","1/h",LOG10},
 { rgN0,"rgN0","1/h",LOG10RATIO},
 {raN0,"raN0","1/h",LOG10RATIO},
 { raNa,"raNa","1/h",LOG10RATIO},
 { rRLN0,"rRLN0","1/h",LOG10RATIO},
 { rRLNa,"rRLNa","1/h",LOG10RATIO},
 { k_NxN,"k_NxN","1/h",LOG10},
 { k_NxAb,"k_NxAb","1/h",LOG10},
 { k_N0NA_A,"k_N0NA_A","1/h",LOG10},
 { KA_N0_Na,"KA_N0_Na","1/h",LOG10},
 { Kg_N0_Na,"Kg_N0_Na","1/h",LOG10},
 {Ka_N0_Na,"Ka_N0_Na","1/h",LOG10},
 { gT0,"gT0","1/h",LOG10},
 { gTs,"gTs","1/h",LOG10},
 { gIT,"gIT","1/h",LOG10},
 { aT0,"aT0","1/h",LOG10},
 { aTs,"aTs","1/h",LOG10},
 { aIT,"aIT","1/h",LOG10},
 { muT0,"muT0","1/h",LOG10},
 { muTs,"muTs","1/h",LOG10},
 { muIT,"muIT","1/h",LOG10},
 { muaT,"muaT","1/h",LOG10},
 {KamuT,"KamuT","1/h",LOG10},
 { kappaT0,"kappaT0","1/h",LOG10},
 { kappaTs,"kappaTs","1/h",LOG10},
 { kappaIT,"kappaIT","1/h",LOG10},
 { rgT0,"rgT0","1/h",LOG10RATIO},
 { rgTs,"rgTs","1/h",LOG10RATIO},
 { raT0,"raT0","1/h",LOG10RATIO},
 { raTs,"raTs","1/h",LOG10RATIO},
 { rRT0,"rRT0","1/h",LOG10RATIO},
 {K_TxAb,"K_TxAb","1/h",LOG10},
 { tAp,"tAp","1/h",LOG10},
 { mua,"mua","1/h",LOG10},
 { mug,"mug","1/h",LOG10},
 { muAg,"muAg","1/h",LOG10},
 { sTym,"sTym","1/h",LOG10},
 { maxCells,"maxCells","1/h",LOG10},
 { k_axAb,"k_axAb","1/h",LOG10},
 {k_gxAb,"k_gxAb","1/h",LOG10}
};

}
std::vector<VariableInformation>  BCM_CD137::externalStateVariables()const
{
  return { {nC0,"nC0","units",LOG10},
 {Ag0,"Ag0","units",LOG10},
 {Ab_R,"Ab_R","units",LOG10},
 {Ab_a,"Ab_a","units",LOG10},
 {Ab_g,"Ab_g","units",LOG10},
 {a_rec,"a_rec","units",LOG10},
 {g_rec,"g_rec","units",LOG10}
};

}
std::vector<VariableInformation>  BCM_CD137::observedStateVariables()const{
  if (false)
  return {{a_obs,"a_obs","units",LINEAR},
 {g_obs,"g_obs","units",LINEAR},
 {rRLA_obs,"rRLA_obs","units",LINEAR},
 {rRLN_obs,"rRLN_obs","units",LINEAR},
 {rRT_obs,"rRT_obs","units",LINEAR},
 {rg_A_obs,"rg_A_obs","units",LINEAR},
 {ra_A_obs,"ra_A_obs","units",LINEAR},
 {rg_N_obs,"rg_N_obs","units",LINEAR},
 {ra_N_obs,"ra_N_obs","units",LINEAR},
 {rg_T_obs,"rg_T_obs","units",LINEAR},
 {ra_T_obs,"ra_T_obs","units",LINEAR},
 {prolifRate_obs,"prolifRate_obs","units",LINEAR},
 {apopRateT_obs,"apopRateT_obs","units",LINEAR},
 {nCells_obs,"nCells_obs","units",LINEAR}
};
else
return {{a_obs,"a_obs","units",LOG10},
{g_obs,"g_obs","units",LOG10},
{rRLA_obs,"rRLA_obs","units",LOG10PERCENT},
{rRLN_obs,"rRLN_obs","units",LOG10PERCENT},
{rRT_obs,"rRT_obs","units",LOG10PERCENT},
{rg_A_obs,"rg_A_obs","units",LOG10PERCENT},
{ra_A_obs,"ra_A_obs","units",LOG10PERCENT},
{rg_N_obs,"rg_N_obs","units",LOG10PERCENT},
{ra_N_obs,"ra_N_obs","units",LOG10PERCENT},
{rg_T_obs,"rg_T_obs","units",LOG10PERCENT},
{ra_T_obs,"ra_T_obs","units",LOG10PERCENT},
{prolifRate_obs,"prolifRate_obs","units",LOG10},
{apopRateT_obs,"apopRateT_obs","units",LOG10PERCENT},
{nCells_obs,"nCells_obs","units",LOG10}
};


}

std::map<int,Transformation*>  BCM_CD137::setTransformations()const{
  return {{LINEAR,new LinearTranformation()},
{LOG10, new Log10Tranformation()},
{LOG10RATIO,new Log10RatioTranformation()},
{LOG10PERCENT,new Log10PercentTranformation()}
};

}

std::map<int,RangeFormat*>  BCM_CD137::setRangeFormat()const{
  return{
  {RANGE_1S,new Range_1S()},
{MEAN_SD,new Mean_SD()},
{MEAN_DB,new Mean_DB()}
};

}



BCM_CD137::BCM_CD137()
{
  ABC_BCM::init();
  push_backExperiment("media",
  {
                        {nC0,1E6},
                        {Ag0,0.0},
                        {Ab_R,0.0},
                        {Ab_a,0.0},
                        {Ab_g,0.0},
                        {a_rec,0.0},
                        {g_rec,0.0}
                      },
  {
                        /*1*/   {a_obs,16.0,2.402,0.722,MEAN_SD},
                        /*2*/  {a_obs,48.0,0.738,0.446,MEAN_SD},
                        /*3*/  {a_obs,119.0,0.176,0.113,MEAN_SD},
                        /*4*/  {g_obs,16.0,0.034,0.034,MEAN_SD},
                        /*5*/   {g_obs,48.0,0.219,0.149,MEAN_SD},
                        /*6*/   {g_obs,119.0,0.146,0.094,MEAN_SD},//
                        /*7*/   {rRLA_obs,0.0,3.0,0.6,MEAN_SD},//
                        /*8*/   {rRLA_obs,16.0,1.29,0.4,MEAN_SD},//
                        /*9*/   {rRLA_obs,119.0,0.89,0.49,MEAN_SD},//
                        /*10*/  {rRLN_obs,0.0,1.54,0.67,MEAN_SD},
                        /*10b*/ {rRLN_obs,24.0,1.68,0.70,MEAN_SD},//
                        /*10c*/ {rRLN_obs,120.0,2.14,0.74,MEAN_SD},//
                        /*11*/  {rRT_obs,0.0,2.0,1.1,MEAN_SD},//
                        /*12*/  {rRT_obs,16.0,1.6,0.7,MEAN_SD},//
                        /*13*/  {rRT_obs,24.0, 3.1, 1.8,MEAN_SD},//
                        /*14*/  {rRT_obs,119.0,5.1,1.6,MEAN_SD},//
                        /*15*/  {rg_A_obs,16.0,2.83,1.16,MEAN_SD},//
                        /*16*/  {rg_A_obs,119.0,2.73,2.05,MEAN_SD},//
                        /*17*/  {ra_A_obs,16.0,5.93,2.73,MEAN_SD},//
                        /*18*/  {ra_A_obs,119.0,4.12,2.91,MEAN_SD},//
                        /*19*/  {rg_N_obs,24.0,3.55,1.43,MEAN_SD},//
                        /*20*/  {ra_N_obs,24.0,3.06,1.44,MEAN_SD},//
                        /*21*/  {rg_T_obs,119,3.4,1.1,MEAN_SD},//
                        /*22*/  {ra_T_obs,119.0,1.54,1.33,MEAN_SD},//
                        /*23*/  {apopRateT_obs,119.0,16.55,4.10,MEAN_SD},//
                        /*24*/  {prolifRate_obs,119,3516.0,935.0,MEAN_SD},//
                        /*25*/  {nCells_obs,24.0,1.5e6,0.75e5,MEAN_SD},//
                        /*26*/  {nCells_obs,119.0,1.5e6,0.75e5,MEAN_SD},//
                      }
                      );



  push_backExperiment("mtb",
  {
                        {nC0,1E6},
                        {Ag0,10.0},
                        {Ab_R,0.0},
                        {Ab_a,0.0},
                        {Ab_g,0.0},
                        {a_rec,0.0},
                        {g_rec,0.0}
                      },{
                        /*27*/  {a_obs,16.0,50.29,8.82,MEAN_SD},//
                        /*28*/  {a_obs,48.0,42.00,7.737,MEAN_SD},//
                        /*29*/  {a_obs,119.0,28.309,5.613,MEAN_SD},//
                        /*30*/  {g_obs,16.0,4.794,0.89,MEAN_SD},//
                        /*31*/  {g_obs,48.0,12.29,1.873,MEAN_SD},//
                        /*32*/  {g_obs,119.0,28.209,5.584,MEAN_SD},//
                        /*33*/   {rRLA_obs,0.0,3.0,0.6,MEAN_SD},//
                        /*34*/   {rRLA_obs,16.0,20.16,7.16,MEAN_SD},
                        /*35*/   {rRLA_obs,119.0,6.17,3.6,MEAN_SD},
                        /*35b*/  {rRLN_obs,0.0,1.54,0.67,MEAN_SD},
                        /*36*/    {rRLN_obs,24.0, 11.32,1.36,MEAN_SD},
                        /*36b*/  {rRLN_obs,120.0,2.31,0.74,MEAN_SD},//
                        /*37*/   {rRT_obs,0.0, 2.0, 1.1,MEAN_SD},
                        /*37b*/   {rRT_obs,16.0, 2.7, 1.9,MEAN_SD},
                        /*38*/   {rRT_obs,24.0, 4.1, 2.3,MEAN_SD},//
                        /*39*/   {rRT_obs,119.0,31.2, 6.9,MEAN_SD},
                        /*40*/    {rg_A_obs,16.0,7.66,3.20,MEAN_SD},//
                        /*40b*/  {rg_A_obs,119.0,2.83,1.16,MEAN_SD},//
                        /*40c*/  {ra_A_obs,16.0,12.96,2.52,MEAN_SD},
                        /*41*/    {ra_A_obs,119.0,7.05,4.03,MEAN_SD},
                        /*42*/     {rg_N_obs,24.0,27.18,2.32,MEAN_SD},
                        /*43*/     {ra_N_obs,24.0,5.81,0.97,MEAN_SD},
                        /*44*/     {ra_T_obs,119.0,4.27,0.59,MEAN_SD},
                        /*45*/     {rg_T_obs,119.0,9.4,1.5,MEAN_SD},
                        /*46*/     {apopRateT_obs,119.0,27.61,2.57,MEAN_SD},
                        /*47*/     {prolifRate_obs,119.0,14173.0,1240.0,MEAN_SD},
                        /*48*/     {nCells_obs,24.0,1.5e6,0.75e5,MEAN_SD},//
                        /*49*/     {nCells_obs,119.0,1.5e6,0.75e5,MEAN_SD},//
                      });

  push_backExperiment("CD137block",
  {
                        {nC0,1E6},
                        {Ag0,10.0},
                        {Ab_R,10.0},
                        {Ab_a,0.0},
                        {Ab_g,0.0},
                        {a_rec,0.0},
                        {g_rec,0.0}
                      },{

                        /*50*/        {a_obs,16.0,61.503,8.527,MEAN_SD},
                        /*51*/        {a_obs,48.0,54.45,8.102,MEAN_SD},
                        /*52*/        {a_obs,119.0,38.86,6.632,MEAN_SD},
                        /*53*/   {g_obs,16.0,8.2280,1.251,MEAN_SD},//
                        /*54*/   {g_obs,48.0,7.911,1.208,MEAN_SD},
                        /*55*/   {g_obs,119.0,13.091,2.24,MEAN_SD},//
                        /*56*/   {rg_A_obs,16.0,14.11,4.02,MEAN_SD},//
                        /*57*/   {rg_A_obs,119.0,5.23,1.41,MEAN_SD},//
                        /*58*/   {ra_A_obs,16.0,41.7,6.99,MEAN_SD},//
                        /*59*/   {ra_A_obs,119.0,12.03,4.13,MEAN_SD},//
                        /*60*/   {ra_N_obs,24.0,8.66,1.48,MEAN_SD},//
                        /*61*/   {rg_N_obs,24.0,36.11,3.2,MEAN_SD},
                        /*62*/   {rg_T_obs,119.0,3.1,1.21,MEAN_SD},
                        /*63*/   {ra_T_obs,119.0,1.93,0.893,MEAN_SD},//
                        /*64*/   {apopRateT_obs,119.0,43.13,3.86,MEAN_SD},//
                        /*65*/   {prolifRate_obs,119.0,5740.0,397.0,MEAN_SD},
                        /*66*/        {nCells_obs,24.0,1.0e6,5.0e5,MEAN_SD},//
                        /*67*/        {nCells_obs,119.0,1.e6,5.0e5,MEAN_SD},//

                      }
                      );

  push_backExperiment("TNFblock",
  {
                        {nC0,1E6},
                        {Ag0,10.0},
                        {Ab_R,0.0},
                        {Ab_a,1.0},
                        {Ab_g,0.0},
                        {a_rec,0.0},
                        {g_rec,0.0}
                      },   {
                        /*68*/    {rRLA_obs,16.0,13.93,4.25,MEAN_SD},//
                        /*69*/    {rRLA_obs,119.0,2.9,1.16,MEAN_SD},//
                        /*70*/    {rRT_obs,119.0,24.3,5.1,MEAN_SD},//

                      });
  push_backExperiment("TNFrecomb",
  {
                        {nC0,1E6},
                        {Ag0,10.0},
                        {Ab_R,0.0},
                        {Ab_a,0.0},
                        {Ab_g,0.0},
                        {a_rec,1.0},
                        {g_rec,0.0}
                      },   {
                        /*71*/    {rRLA_obs,16.0,26.62,7.32,MEAN_SD},//
                        /*72*/    {rRLA_obs,119.0,11.9,3.53,MEAN_SD},//
                        /*73*/    {rRT_obs,119.0,28.3,6.9,MEAN_SD},//

                      });
  push_backExperiment("IFNblock",
  {
                        {nC0,1E6},
                        {Ag0,10.0},
                        {Ab_R,0.0},
                        {Ab_a,0.0},
                        {Ab_g,1.0},
                        {a_rec,0.0},
                        {g_rec,0.0}
                      },   {
                        /*74*/    {rRLA_obs,16.0,65.17,11.43,MEAN_SD},//
                        /*75*/    {rRLA_obs,119.0,49.84,13.72,MEAN_SD},//
                        /*76*/    {rRT_obs,119.0,36.2,8.4,MEAN_SD},//

                      });

  push_backExperiment("IFNrecomb",
  {
                        {nC0,1E6},
                        {Ag0,10.0},
                        {Ab_R,0.0},
                        {Ab_a,0.0},
                        {Ab_g,0.0},
                        {a_rec,0.0},
                        {g_rec,1.0}
                      },   {
                        /*77*/    {rRLA_obs,16.0,13.77,4.1,MEAN_SD},//
                        /*78*/    {rRLA_obs,119.0,5.59,3.11,MEAN_SD},//
                        /*79*/    {rRT_obs,119.0,32.6,6.6,MEAN_SD},//

                      });






  push_backParameters("PriorsIniciales",
  {
                        {rA,0.7,0.88,RANGE_1S},
                        {rN,1.22E-2,3.55E-2,RANGE_1S},
                        {rT,0.58,0.82,RANGE_1S},
                        {rTsp,2.00E-4,2.00E-3,RANGE_1S},
                        {gA0,8.00E-7,2.50E-6,RANGE_1S},
                        {gAa,8.00E-7,2.50E-6,RANGE_1S},
                        {gIA,0.100,10.0,RANGE_1S},
                        {aA0,6.40E-7,6.40E-4,RANGE_1S},
                        {aAa,6.40E-7,6.40E-4,RANGE_1S},
                        {aIA,0.100,10.0,RANGE_1S},
                        {muA0,2.80E-5, 2.80E-3,RANGE_1S},
                        {muAa,1.50E-4, 1.50E-2,RANGE_1S},
                        {muIA,0.100,10.0,RANGE_1S},
                        {muaA,4.20E-4,4.20E-2,RANGE_1S},
                        {KamuA,1.00E-3,10.0,RANGE_1S},
                        {kappaA,4.00E-12,0.400,RANGE_1S},
                        {rgA0,1.00E-2,0.100,RANGE_1S},
                        {rgAa,0.25,0.75,RANGE_1S},
                        {raA0,1.00E-2,0.100,RANGE_1S},
                        {rRLA0,1.00E-2,5.0E-2,RANGE_1S},
                        {k_AxA,1.00E-6,1.00E-2,RANGE_1S},
                        {k_AxN,1.00E-6,1.00E-2,RANGE_1S},
                        {k_AxT,1.00E-6,1.00E-2,RANGE_1S},
                        {k_AxAg,1.00E-8,1.00E-5,RANGE_1S},
                        {k_AxAg_ag,1.00E-3,0.100,RANGE_1S},
                        {k_AxAb,5.00E-2,5.00,RANGE_1S},
                        {Kg_A0Aa,3.00E-2,60.0,RANGE_1S},
                        {Ka_A0Aa,3.40E-2,150.0,RANGE_1S},
                        {gN0,4.2E-8,4.2E-4,RANGE_1S},
                        {gNa,1.30E-4,0.130,RANGE_1S},
                        {gIN,2.0,10.0,RANGE_1S},
                        {aN0,3.75E-10,3.75E-8,RANGE_1S},
                        {aNa,3.75E-10,3.75E-8,RANGE_1S},
                        {aIN,2.0,10.0,RANGE_1S},
                        {muN0,4.00E-8,4.00E-4,RANGE_1S},
                        {muNa,4.00E-8,4.00E-4,RANGE_1S},
                        {muaN,4.20E-4,4.20E-1,RANGE_1S},
                        {KamuN,3.40E-2,150.0,RANGE_1S},
                        {kappaN0,3.00E-5,3.00E-3,RANGE_1S},
                        {kappaNa,3.00E-5,3.00E-3,RANGE_1S},
                        {rgN0,3.00E-2,5.00E-2,RANGE_1S},
                        {raN0,1.00E-2, 5.00E-2,RANGE_1S},
                        {raNa,0.10,0.50,RANGE_1S},
                        {rRLN0,1.00E-2,0.170,RANGE_1S},
                        {rRLNa,9.00E-2,0.490,RANGE_1S},
                        {k_NxN,1.00E-6,1.00E-2,RANGE_1S},
                        {k_NxAb,5.00E-2,5.00,RANGE_1S},
                        {k_N0NA_A ,1.00E-6,1.00E-3,RANGE_1S},
                        {KA_N0_Na,2.00E3,2.00E5,RANGE_1S},
                        {Kg_N0_Na,3.00E-2,60.0,RANGE_1S},// de mas
                        {Ka_N0_Na,3.40E-2,150.0,RANGE_1S}, //de mas
                        {gT0,8.33E-7,2.75E-6,RANGE_1S},
                        {gTs,1.20E-6,1.20E-3,RANGE_1S},
                        {gIT,1.20,6.70,RANGE_1S},
                        {aT0,2.50E-11,1.00E-8,RANGE_1S},
                        {aTs,2.50E-10,1.00E-7,RANGE_1S},
                        {aIT,1.20,3.50,RANGE_1S},
                        {muT0,1.00E-2,5.20E-2,RANGE_1S},
                        {muTs,8.30E-3,8.30E-2,RANGE_1S},
                        {muIT,1.00,5.00,RANGE_1S},
                        {muaT,4.20E-4,4.20E-2,RANGE_1S},
                        {KamuT,3.40E-2,150.0,RANGE_1S},
                        {kappaT0,1.66E-4,1.66E-2,RANGE_1S},
                        {kappaTs,0.0083,0.83,RANGE_1S},
                        {kappaIT,1.00,12.0,RANGE_1S},
                        {rgT0,1.00E-2,2.00E-2,RANGE_1S},
                        {rgTs,0.396,0.520,RANGE_1S},
                        {raT0,1.00E-2,5.00E-2,RANGE_1S},
                        {raTs,0.201,0.278,RANGE_1S},
                        {rRT0,5.00E-2,0.200,RANGE_1S},
                        {K_TxAb,0.100,2.00,RANGE_1S},
                        {tAp,2.00,20.0,RANGE_1S},
                        {mua,4.10E-2,0.150,RANGE_1S},
                        {mug,4.10E-2,0.150,RANGE_1S},
                        {muAg,4.10E-2,0.150,RANGE_1S},
                        {sTym,1.00E-4,0.100,RANGE_1S},
                        {maxCells,5.00E5,3.00E6,RANGE_1S},
                        {k_axAb,5.00E-2,5.00,RANGE_1S},
                        {k_gxAb,5.00E-2,5.00,RANGE_1S},

                      }
                      );

  setPrior("PriorsIniciales");


}



std::vector<double> BCM_CD137::internalStatesStart(const std::vector<double>& p,
                                                   const std::vector<double>& e)const
{
  std::vector<double>i(interNames.size());

  i[A_0]=(1.0-p[rT])*p[rA]*e[nC0];
  i[A_a]=0.0;
  i[A_s]=0.0;
  i[A_Ab]=0.0;
  i[A_s_Ab]=0.0;

  i[N_0]=(1.0-p[rT])*(1.0-p[rA])*e[nC0];
  i[N_a]=0.0;
  i[N_s]=0.0;
  i[N_Ab]=0.0;
  i[N_s_Ab]=0.0;

  i[T_ns]=p[rT]*(1.0-p[rTsp])*e[nC0];
  i[T_0]=p[rT]*(p[rTsp])*e[nC0];
  i[T_s]=0.0;
  i[T_bl]=0.0;

  i[a]=e[a_rec];
  i[g]=e[g_rec];
  i[Ag]=e[Ag0];


  return i;
}

std::vector<double> BCM_CD137::internalStatesDerivative(const std::vector<double>& p,
                                                        const std::vector<double>& e,
                                                        const std::vector<double>& i)const
{
  std::vector<double> D(interNames.size());

  double nC=i[A_0]+i[ A_a ]+i[ A_s ]+i[ A_Ab ]+i[ A_s_Ab  ]+
      i[N_0]+i[ N_a]+i[N_s ]+i[N_Ab]+i[N_s_Ab]+
      i[T_ns]+i[T_0]+i[T_s]+i[T_bl];

  double r_kappa=(p[maxCells]-nC)/p[maxCells];
  D[A_0]=-p[muA0]*i[A_0]
      -p[k_AxAg]*i[Ag]*i[A_0]
      -p[k_AxAg_ag]*(i[g]/(i[g]+p[Kg_A0Aa]))*(i[a]/(i[a]+p[Ka_A0Aa]))*i[A_0];




  D[A_a]=p[k_AxAg]*i[Ag]*i[A_0]
      +p[k_AxAg_ag]*(i[g]/(i[g]+p[Kg_A0Aa]))*(i[a]/(i[a]+p[Ka_A0Aa]))*i[A_0]
      -p[muAa]*i[A_a]
      -p[muaA]*(i[a]/(i[a]+p[KamuA]))*i[A_a]
      -p[k_AxA]*(2.0*i[A_a]+i[A_s])*i[A_a]
      -p[k_AxN]*p[rRLNa]*(i[N_a]+i[N_s])*i[A_a]
      -p[k_AxT]*p[K_TxAb]/(p[K_TxAb]+e[Ab_R])*i[T_0]*i[A_a]
      -p[k_AxAb]*e[Ab_R]*i[A_a];


  D[A_s]=p[k_AxA]*(2.0*i[A_a]+i[A_s])*i[A_a]
      +p[k_AxN]*p[rRLNa]*(i[N_a]+i[N_s])*i[A_a]
      +p[k_AxT]*p[K_TxAb]/(p[K_TxAb]+e[Ab_R])*i[T_0]*i[A_a]
      -p[k_AxAb]*e[Ab_R]*i[A_s]
      +p[kappaA]*r_kappa*i[A_s]
      -p[muAa]*p[muIA]*i[A_s]
      -p[muaA]*(i[a]/(i[a]+p[KamuA]))*i[A_s];




  D[A_Ab]=p[k_AxAb]*e[Ab_R]*i[A_a]
      -p[muAa]*i[A_Ab]
      -p[muaA]*(i[a]/(i[a]+p[KamuA]))*i[A_Ab]
      -p[k_AxA]*(2.0*i[A_a]+i[A_s])*i[A_Ab]
      -p[k_AxN]*p[rRLNa]*(i[N_a]+i[N_s])*i[A_Ab]
      -p[k_AxT]*p[K_TxAb]/(p[K_TxAb]+e[Ab_R])*i[T_0]*i[A_Ab];




  D[A_s_Ab]=p[k_AxAb]*e[Ab_R]*i[A_s]
      +p[k_AxA]*(2.0*i[A_a]+i[A_s])*i[A_Ab]
      +p[k_AxN]*p[rRLNa]*(i[N_a]+i[N_s])*i[A_Ab]
      +p[k_AxT]*p[K_TxAb]/(p[K_TxAb]+e[Ab_R])*i[T_0]*i[A_Ab]
      -p[muAa]*p[muIA]*i[A_s_Ab]
      -p[muaA]*(i[a]/(i[a]+p[KamuA]))*i[A_s_Ab]
      +p[kappaA]*r_kappa*i[A_s_Ab];




  D[N_0]=p[kappaN0]*r_kappa*i[N_0]
      -p[muN0]*i[N_0]
      -(p[k_N0NA_A]*(i[A_a]+i[A_Ab]+p[aIA]*(i[A_s]+i[A_s_Ab]))
        /(i[A_a]+i[A_Ab]+p[aIA]*(i[A_s]+i[A_s_Ab])+p[KA_N0_Na])
        *(i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab])*i[Ag])*i[N_0];


  //  double n0_a_NA=(p[k_N0NA_A]*(i[A_a]+i[A_Ab]+p[aIA]*(i[A_s]+i[A_s_Ab]))
  //                  /(i[A_a]+i[A_Ab]+p[aIA]*(i[A_s]+i[A_s_Ab])+p[KA_N0_Na])
  //                  *(i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab])*i[Ag])*i[N_0];


  D[N_a]=p[kappaNa]*r_kappa*i[N_a]
      +(p[k_N0NA_A]*(i[A_a]+i[A_Ab]+p[aIA]*(i[A_s]+i[A_s_Ab]))
        /(i[A_a]+i[A_Ab]+p[aIA]*(i[A_s]+i[A_s_Ab])+p[KA_N0_Na])
        *(i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab])*i[Ag])*i[N_0]
      -p[k_NxN]*p[rRLNa]*(p[rRLNa]*2.0*i[N_a]+i[N_Ab]+i[N_s]+i[N_s_Ab])*i[N_a]
      -p[k_AxN]*p[rRLNa]*(i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab])*i[N_a]
      -p[k_NxAb]*e[Ab_R]*i[N_a]
      -p[muNa]*i[N_a]
      -p[muaN]*(i[a]/(i[a]+KamuN))*i[N_a];


  //double Na_aNs=p[k_NxN]*p[rRLNa]*(p[rRLNa]*2.0*i[N_a]+i[N_Ab]+i[N_s]+i[N_s_Ab])*i[N_a]
  //    +p[k_AxN]*p[rRLNa]*(i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab])*i[N_a];


  D[N_s]=p[kappaNa]*r_kappa*i[N_s]
      +p[k_NxN]*p[rRLNa]*(p[rRLNa]*2.0*i[N_a]+i[N_Ab]+i[N_s]+i[N_s_Ab])*i[N_a]
      +p[k_AxN]*p[rRLNa]*(i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab])*i[N_a]
      -p[muNa]*i[N_s]
      -p[muaN]*(i[a]/(i[a]+p[KamuN]))*i[N_s]
      -p[k_NxAb]*e[Ab_R]*i[N_s];





  D[N_Ab]=p[kappaNa]*r_kappa*i[N_Ab]
      +p[k_NxAb]*e[Ab_R]*i[N_s]
      -p[muNa]*i[N_Ab]
      -p[muaN]*(i[a]/(i[a]+p[KamuN]))*i[N_Ab];




  D[N_s_Ab]=p[kappaNa]*r_kappa*i[N_s_Ab]
      +p[k_NxAb]*e[Ab_R]*i[N_s]
      -p[muNa]*i[N_s_Ab]
      -p[muaN]*(i[a]/(i[a]+p[KamuN]))*i[N_s_Ab];




  D[T_ns]=p[kappaT0]*r_kappa*i[T_ns]
      -p[muT0]*i[T_ns];





  D[T_0]=p[kappaT0]*r_kappa*i[T_0]
      -p[muT0]*i[T_0]
      -p[k_AxT]*(i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab])*i[T_0];




  D[T_s]=p[kappaTs]*r_kappa*i[T_s]
      -p[muTs]*i[T_s]
      -p[muaT]*i[a]/(i[a]+p[KamuT])*i[T_s]
      +p[k_AxT]*p[K_TxAb]/(p[K_TxAb]+e[Ab_R])*(i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab])*i[T_0];


  D[T_bl]=p[kappaTs]*p[kappaIT]*r_kappa*i[T_bl]
      -p[muTs]/p[muIT]*i[T_bl]
      -p[muaT]*i[a]/(i[a]+p[KamuT])*i[T_bl]
      +p[k_AxT]*e[Ab_R]/(p[K_TxAb]+e[Ab_R])*(i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab])*i[T_0];





  D[a]=p[aA0]*p[raA0]*i[A_0]
      +p[aAa]*(i[A_a]+i[A_Ab])
      +p[aAa]*p[aIA]*(i[A_s]+i[A_s_Ab])
      +p[aN0]*p[raN0]*i[N_0]
      +p[aNa]*p[raNa]*(i[N_a]+i[N_Ab])
      +p[aNa]*p[raNa]/p[aIN]*(i[N_s]+i[N_s_Ab])
      +p[aT0]*p[raT0]*(i[T_0]+i[T_ns])
      +p[aTs]*p[raTs]*i[T_s]
      +p[aTs]*p[raTs]*p[aIT]*i[T_bl]
      -p[mua]*i[a]
      -p[k_axAb]*i[a]*e[Ab_a];


  D[g]=p[gA0]*p[rgA0]*i[A_0]
      +p[gAa]*p[rgAa]*(i[A_a]+i[A_Ab])
      +p[gAa]*p[gIA]*p[rgAa]*(i[A_s]+i[A_s_Ab])
      +p[gN0]*p[rgN0]*i[N_0]
      +p[gNa]*(i[N_a]+i[N_Ab])
      +p[gNa]/p[gIN]*(i[N_s]+i[N_s_Ab])
      +p[gT0]*p[rgT0]*(i[T_0]+i[T_ns])
      +p[gTs]*p[rgTs]*i[T_s]
      +p[gTs]*p[rgTs]*p[gIT]*i[T_bl]
      -p[mug]*i[g]
      -p[k_gxAb]*i[g]*e[Ab_g];

  D[Ag]=-p[muAg]*i[Ag];


  return D;
}
double BCM_CD137::observedVariablesEquations(const std::vector<double>& p,
                                             const std::vector<double>& e,
                                             const std::vector<double>& i,
                                             int observedIndex,
                                             double currentsum,
                                             double dt,
                                             bool finalcalculation)const
{
  switch (observedIndex)
    {  case a_obs:
      return i[a];
    case g_obs:
      return i[g];
    case rRLA_obs:
      {
        double A_tot=i[A_0]+i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab];
        return (p[rRLA0]*i[A_0]+i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab])*100.0/A_tot;
      }
    case rRLN_obs:
      {
        double N_tot=i[N_0]+i[N_a]+i[N_Ab]+i[N_s]+i[N_s_Ab];
        return (p[rRLN0]*i[N_0]+p[rRLNa]*(i[N_a]+i[N_Ab]+i[N_s]+i[N_s_Ab]))*100.0/N_tot;
      }
    case rRT_obs:
      {
        double T_tot=i[T_ns]+i[T_0]+i[T_s]+i[T_bl];
        return (p[raT0]*(i[T_ns]+i[T_0])+i[T_s]+i[T_bl])*100.0/T_tot;
      }
    case rg_A_obs:
      {
        double A_tot=i[A_0]+i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab];
        return (p[rgA0]*i[A_0]+p[rgAa]*(i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab]))*100.0/A_tot;
      }
    case ra_A_obs:
      {
        double A_tot=i[A_0]+i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab];
        return (p[raA0]*i[A_0]+i[A_a]+i[A_Ab]+i[A_s]+i[A_s_Ab])*100.0/A_tot;
      }
    case rg_N_obs:
      {
        double N_tot=i[N_0]+i[N_a]+i[N_Ab]+i[N_s]+i[N_s_Ab];
        return (p[rgN0]*i[N_0]+i[N_a]+i[N_Ab]+i[N_s]+i[N_s_Ab])*100.0/N_tot;
      }
    case ra_N_obs:
      {
        double N_tot=i[N_0]+i[N_a]+i[N_Ab]+i[N_s]+i[N_s_Ab];
        return (p[raN0]*i[N_0]+p[raNa]*(i[N_a]+i[N_Ab]+i[N_s]+i[N_s_Ab]))*100.0/N_tot;
      }
    case rg_T_obs:
      {
        double T_tot=i[T_ns]+i[T_0]+i[T_s]+i[T_bl];
        return (p[rgT0]*(i[T_ns]+i[T_0])+p[rgTs]*(i[T_s]+i[T_bl]))*100.0/T_tot;
      }

    case ra_T_obs:
      {
        double T_tot=i[T_ns]+i[T_0]+i[T_s]+i[T_bl];
        return (p[raT0]*(i[T_ns]+i[T_0])+p[raTs]*(i[T_s]+i[T_bl]))*100.0/T_tot;
      }
    case prolifRate_obs:
      {
        double nC=i[A_0 ]+i[ A_a ]+i[ A_s ]+i[ A_Ab ]+i[ A_s_Ab  ]+
            i[N_0]+i[ N_a]+i[N_s ]+i[N_Ab]+i[N_s_Ab]+
            i[T_ns]+i[T_0]+i[T_s]+i[T_bl];

        double r_kappa=(p[maxCells]-nC)/p[maxCells];

        if (!finalcalculation)
          return currentsum+dt*(p[kappaA]*(i[A_s]+i[A_s_Ab])+
                                p[kappaN0]*i[N_0]+p[kappaNa]*(i[N_a]+i[N_s]+i[N_Ab]+i[N_s_Ab])+
                                p[kappaT0]*(i[T_ns]+i[T_0])+
                                p[kappaTs]*(i[T_s]+p[kappaIT]*i[T_bl]))*r_kappa;
        else

          return p[sTym]*(currentsum+dt*(p[kappaA]*(i[A_s]+i[A_s_Ab])+
                                         p[kappaN0]*i[N_0]+p[kappaNa]*(i[N_a]+i[N_s]+i[N_Ab]+i[N_s_Ab])+
                                         p[kappaT0]*(i[T_ns]+i[T_0])+
                                         p[kappaTs]*(i[T_s]+p[kappaIT]*i[T_bl]))*r_kappa);
      }
    case apopRateT_obs:
      {
        if (!finalcalculation)
          {
            return currentsum+dt*(p[muT0]*(i[T_ns]+i[T_0])+
                                  p[muTs]*(i[T_s]+p[muIT]*i[T_bl])+
                                  p[muaT]*i[a]/(i[a]+p[KamuT])*(i[T_s]+i[T_bl]));
          }
        else
          {
            double T_tot=i[T_ns]+i[T_0]+i[T_s]+i[T_bl];
            double Apt=(currentsum+
                        dt*(p[muT0]*(i[T_ns]+i[T_0])+
                            p[muTs]*(i[T_s]+p[muIT]*i[T_bl])+
                            p[muaT]*i[a]/(i[a]+p[KamuT])*(i[T_s]+i[T_bl])));

            return 100.0*Apt/(T_tot+Apt);
          }
      }
    case nCells_obs:
      return i[A_0 ]+i[ A_a ]+i[ A_s ]+i[ A_Ab ]+i[ A_s_Ab  ]+
          i[N_0]+i[ N_a]+i[N_s ]+i[N_Ab]+i[N_s_Ab]+
          i[T_ns]+i[T_0]+i[T_s]+i[T_bl];


    }
}



std::string BCM_CD137::identifierName()const
{
  return "BCM_CD137";
}


double BCM_CD137::observedVariableIntegrationTime(const std::vector<double>& p,
                                                  const std::vector<double>& e,
                                                  int observedIndex)const
{
  switch (observedIndex) {
    case prolifRate_obs:
      return 16.0;
    case apopRateT_obs:
      return p[tAp];
    default:
      return 0.0;
    }
}



