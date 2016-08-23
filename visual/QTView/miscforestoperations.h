#ifndef MISCFORESTOPERATIONS_H
#define MISCFORESTOPERATIONS_H
#include "parser.h"
#include <QLabel>
#include <QGraphicsScene>

//class MiscForestOperations
//{
//public:
//    MiscForestOperations();
//};

QVector<int> *setupRelation(QStringList &list);

QVector<int> *setupSet(QStringList &list);

QVector<int> *setupForest1(Parser *parser, int &f1base
                           , int &f1PenHeight
                           , QString &f1Name);

QVector<int> *setupForest2(Parser *parser
                           , int &f2base
                           , int &f2PenHeight
                           , QString &f2Name);

void playForests(Parser *parser
                 , QLabel *labelForest1
                 , QLabel *labelForest2
                 , QVector<int> *&forest1
                 , QVector<int> *&forest2
                 , QGraphicsScene *&forest1Scene
                 , QGraphicsScene *&forest2Scene
                 , QVector<int> &isInForest1Updates
                 , QVector<int> &isInForest2Updates
                 , QVector<int> &Forest1Updates
                 , QVector<int> &Forest2Updates
                 , int f1base, int f2base
                 , int f1PenHeight, int f2PenHeight
                 , int &oldSeconds
                 , int &oldMiliseconds
                 , QString f1Name, QString f2Name
                 , bool &endOfFile
                 , int f1HorizontalReductionValue
                 , int f2HorizontalReductionValue
                 , int speedFactorValue);

void redrawForest(QGraphicsScene *forestScene
                   , QVector<int> *forest, int base
                  , int penHeight, int reductionValue);

//void redrawForest2(QGraphicsScene *forest2Scene);
#endif // MISCFORESTOPERATIONS_H
