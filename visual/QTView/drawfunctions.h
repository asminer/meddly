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

void drawUpdates(QVector<int> *&forest
                 , QVector<int> &isInForestUpdates
                 , QVector<int> &ForestUpdates
                 , QGraphicsScene *&forestScene
                 , int fPenHeight, int fbase, int fHorizontalReductionValue
                 );

void drawRulers(QGraphicsScene *forestScene
                , int penHeight, int forestsize, int setOrRelation);

#endif // DRAWFUNCTIONS_H
