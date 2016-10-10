#ifndef DRAWFUNCTIONS_H
#define DRAWFUNCTIONS_H

#include <QGraphicsScene>
#include <QVector>
//class Drawfunctions
//{
//public:
//    Drawfunctions();
//};

void drawSet(QGraphicsScene *forestScene
             , QVector<int> *forest, int penHeight
             , int reductionValue);

void drawRelation(QGraphicsScene *forestScene
                  , QVector<int> *forest, int penHeight
                  , int reductionValue);

int generatePenHeight(int size);

void drawUpdates(QVector<int> *&forest1
                 , QVector<int> *&forest2
                 , QVector<int> &isInForest1Updates
                 , QVector<int> &isInForest2Updates
                 , QVector<int> &Forest1Updates
                 , QVector<int> &Forest2Updates
                 , QGraphicsScene *&forest1Scene
                 , QGraphicsScene *&forest2Scene
                 , int f1PenHeight, int f2PenHeight
                 , int f1base, int f2base
                 , int f1HorizontalReductionValue
                 , int f2HorizontalReductionValue);

void drawRulers(QGraphicsScene *forestScene
                , int penHeight, int forestsize, int setOrRelation);

#endif // DRAWFUNCTIONS_H
