#include "drawfunctions.h"
#include "timefunctions.h"

#include <QMessageBox>
#include <QGraphicsLineItem>
#include <QGraphicsItem>

void drawSet(QGraphicsScene *forestScene
             , QVector<int> *forest
             , int penHeight, int reductionValue)
{
    QPen green(Qt::green);
    green.setWidth(penHeight);
    int startLine = 0;
    int endLine = 0;

    int g = forest->size();
    int y = 0;

    for (int var = 0; var < g; ++var)
    {
        endLine = ( (forest->at(var))/2 ) /reductionValue;
        startLine = 0 - endLine;

        y = var *penHeight;

//        QGraphicsLineItem *line = new QGraphicsLineItem(startLine
//                                                       , y
//                                                       , endLine
//                                                       , y);
//        line->setPen(green);
//        forestScene->addItem(line);

        forestScene->addLine(startLine
                             , y
                             , endLine
                             , y
                             , green);

    }
}

void drawRelation(QGraphicsScene *forestScene
                  , QVector<int> *forest
                  , int penHeight, int reductionValue)
{
    QPen green(Qt::green);
    green.setWidth(penHeight);

    QPen darkGreen(Qt::darkGreen);
    darkGreen.setWidth(penHeight);
    int startLine = 0;
    int endLine = 0;

    int posYlevel = 0;
    int negYlevel = 0;

    int top = forest->size() - 1;
    int middle = top/2;

    for (int var = 0; var < middle; ++var)
    {
        //Positive
        endLine = ( (forest->at(middle + var + 1))/2 ) / reductionValue;
        startLine = 0 - endLine;

        posYlevel = ((var*2)+2) *penHeight;

        forestScene->addLine(startLine
                             , posYlevel
                             , endLine
                             , posYlevel
                             , green);

        //Negitive
        endLine = ( (forest->at(middle - var -1))/2 ) / reductionValue;
        startLine = 0 - endLine;
        negYlevel = ((var*2)+1) * penHeight;

        forestScene->addLine(startLine
                             , negYlevel
                             , endLine
                             , negYlevel
                             , darkGreen);
        delayMili(1);
    }
}

int generatePenHeight(int size)
{
    if(size <20)
    {
        return 8;
    }else if(size >= 20 && size <= 100)
    {
        return 5;
    }else
    {
        return 1;
    }
}

void drawUpdates(QVector<int> *&forest1
                 , QVector<int> *&forest2
                 , QVector<int> &isInForest1Updates
                 , QVector<int> &isInForest2Updates
                 , QVector<int> &Forest1Updates
                 , QVector<int> &Forest2Updates
                 , QGraphicsScene *&forest1Scene
                 , QGraphicsScene *&forest2Scene
                 , int f1PenHeight
                 , int f2PenHeight
                 , int f1base
                 , int f2base
                 , int f1HorizontalReductionValue
                 , int f2HorizontalReductionValue)
{
    QPen green(Qt::green);
    QPen darkgreen(Qt::darkGreen);

    QGraphicsItem *line = NULL;
    QGraphicsLineItem *lineI = NULL;
    qreal x = 0;
    qreal y = 0;
    if(f1base > 0)
    {
        while (!Forest1Updates.isEmpty())
        {
            int indexOfForest = Forest1Updates.takeFirst();
            isInForest1Updates.replace(indexOfForest, 0);
            y = indexOfForest * f1PenHeight;
            line = forest1Scene->itemAt(x, y, QTransform());
            lineI = qgraphicsitem_cast<QGraphicsLineItem*>(line);
            //forest1Scene->removeItem(line);
            //delete line;

            int endline = forest1->at(indexOfForest)/2;
            endline = endline/f1HorizontalReductionValue;
            int startline = 0 - endline;

            if(lineI == NULL)
            {
                green.setWidth(f1PenHeight);
                forest1Scene->addLine(startline
                                      , y
                                      , endline
                                      , y
                                      , green
                                      );
            }else
            {
                lineI->setLine(startline
                               , y
                               , endline
                               , y);
            }
            //delayMili(1);
        }
    }else
    {
        while (!Forest1Updates.isEmpty())
        {
            int indexOfForest = Forest1Updates.takeFirst();
            isInForest1Updates.replace(indexOfForest, 0);
            //positive levels
            int middle = forest1->size()/2;
            if(indexOfForest > middle)
            {
                int indexOnScene = (indexOfForest - middle) * 2;

                y = indexOnScene * f1PenHeight;
                line = forest1Scene->itemAt(x, y, QTransform());
                forest1Scene->removeItem(line);
                delete line;

                int endline = forest1->at(indexOfForest)/2;
                endline = endline/f1HorizontalReductionValue;
                int startline = 0 - endline;

                green.setWidth(f1PenHeight);
                forest1Scene->addLine(startline
                                      , y
                                      , endline
                                      , y
                                      , green
                                      );
                //delayMili(1);

            }else//negitive levels
            {
                int indexOnScene = (forest1->size() - 1)
                        - (indexOfForest * 2) -1;
                y = indexOnScene * f1PenHeight;
                line = forest1Scene->itemAt(x, y, QTransform());
                forest1Scene->removeItem(line);
                delete line;

                int endline = forest1->at(indexOfForest)/2;
                endline = endline/f1HorizontalReductionValue;
                int startline = 0 - endline;

                darkgreen.setWidth(f1PenHeight);
                forest1Scene->addLine(startline
                                      , y
                                      , endline
                                      , y
                                      , darkgreen
                                      );
                //delayMili(1);
            }
        }
    }


    if(f2base > 0)
    {
        while (!Forest2Updates.isEmpty())
        {
            int indexOfForest = Forest2Updates.takeFirst();
            isInForest2Updates.replace(indexOfForest, 0);
            y = indexOfForest * f2PenHeight;
            line = forest2Scene->itemAt(x, y, QTransform());
            forest2Scene->removeItem(line);
            delete line;

            int endline = forest2->at(indexOfForest)/2;
            endline = endline/f2HorizontalReductionValue;
            int startline = 0 - endline;

            green.setWidth(f2PenHeight);
            forest2Scene->addLine(startline
                                  , y
                                  , endline
                                  , y
                                  , green
                                  );
            //delayMili(1);
        }

    }else
    {
        while (!Forest2Updates.isEmpty())
        {
            int indexOfForest = Forest2Updates.takeFirst();
            isInForest2Updates.replace(indexOfForest, 0);

            //positive levels
            int middle = forest2->size()/2;
            if(indexOfForest > middle)
            {
                int indexOnScene = (indexOfForest - middle) * 2;

                y = indexOnScene * f2PenHeight;
                line = forest2Scene->itemAt(x, y, QTransform());
                forest2Scene->removeItem(line);
                delete line;

                int endline = forest2->at(indexOfForest)/2;
                endline = endline/f2HorizontalReductionValue;
                int startline = 0 - endline;

                green.setWidth(f2PenHeight);
                forest2Scene->addLine(startline
                                      , y
                                      , endline
                                      , y
                                      , green
                                      );
                //delayMili(1);

            }else//negitive levels
            {
                int indexOnScene = (forest2->size() - 1)
                        - (indexOfForest * 2) -1;
                y = indexOnScene * f2PenHeight;
                line = forest2Scene->itemAt(x, y, QTransform());
                forest2Scene->removeItem(line);
                delete line;

                int endline = forest2->at(indexOfForest)/2;
                endline = endline/f2HorizontalReductionValue;
                int startline = 0 - endline;

                darkgreen.setWidth(f2PenHeight);
                forest2Scene->addLine(startline
                                      , y
                                      , endline
                                      , y
                                      , darkgreen
                                      );
                //delayMili(1);
            }
        }
    }
}
