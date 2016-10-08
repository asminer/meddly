#include "drawfunctions.h"
#include "timefunctions.h"

#include <QMessageBox>
#include <QGraphicsLineItem>
#include <QGraphicsItem>

/*
 * This functions will add all lines to the scene for sets.
 * This will go through the entire forest and add the lines left aligned, starting on the 0 axis.
 *
 */
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
        //endLine = ( (forest->at(var))/2 ) /reductionValue;
        //startLine = 0 - endLine;

        endLine = forest->at(var) /reductionValue;

        y = var *penHeight;

        forestScene->addLine(0//startLine
                             , y
                             , endLine
                             , y
                             , green);

    }
}

/*
 * This functions will add all lines to the scene for relations.
 * This will go through the entire forest and add the lines left aligned, starting on the 0 axis.
 *
 */
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
        //Positive lines
        endLine = (forest->at(middle + var + 1)) / reductionValue;
        //startLine = 0 - endLine;

        posYlevel = ((var*2)+2) *penHeight;

        forestScene->addLine(0//startLine
                             , posYlevel
                             , endLine
                             , posYlevel
                             , green);

        //Negitive lines
        //endLine = ( (forest->at(middle - var -1))/2 ) / reductionValue;
        //startLine = 0 - endLine;

        endLine = forest->at(middle - var -1) / reductionValue;

        negYlevel = ((var*2)+1) * penHeight;

        forestScene->addLine(0 //startLine
                             , negYlevel
                             , endLine
                             , negYlevel
                             , darkGreen);
        delayMili(1);
    }
}

/*
 * This functions chooses the default pen height for use later when adding lines to a scene.
 */
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

/*
 * This will update the forest scenes.
 *
 */
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

            int endline = forest1->at(indexOfForest);
            endline = endline/f1HorizontalReductionValue;
//          int startline = 0 - endline;

            //If you try to grab a 0 to 0 line then the
            //pointer will be null. This checks to see if
            //you tried to pull an empty line from the
            //scene. If you did then it will create a new line
            //to add to the appropriate spot.
            if(lineI == NULL)
            {
                green.setWidth(f1PenHeight);
                forest1Scene->addLine(0//startline
                                      , y
                                      , endline
                                      , y
                                      , green
                                      );
            }else//This changes the positionof the
                //end of the line object.
            {
                lineI->setLine(0 //startline
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

                int endline = forest1->at(indexOfForest);
                endline = endline/f1HorizontalReductionValue;
                //int startline = 0 - endline;

                green.setWidth(f1PenHeight);
                forest1Scene->addLine(0//startline
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

                int endline = forest1->at(indexOfForest);
                endline = endline/f1HorizontalReductionValue;
                //int startline = 0 - endline;

                darkgreen.setWidth(f1PenHeight);
                forest1Scene->addLine(0 //startline
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

            int endline = forest2->at(indexOfForest);
            endline = endline/f2HorizontalReductionValue;
            //int startline = 0 - endline;

            green.setWidth(f2PenHeight);
            forest2Scene->addLine(0 //startline
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

                int endline = forest2->at(indexOfForest);
                endline = endline/f2HorizontalReductionValue;
                //int startline = 0 - endline;

                green.setWidth(f2PenHeight);
                forest2Scene->addLine(0 //startline
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

                int endline = forest2->at(indexOfForest);
                endline = endline/f2HorizontalReductionValue;
                //int startline = 0 - endline;

                darkgreen.setWidth(f2PenHeight);
                forest2Scene->addLine(0 //startline
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
