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
    //int startLine = 0;
    int endLine = 0;

    int top = forest->size();
    int y = 0;

    drawRulers(forestScene, penHeight
               , top, 0);

    for (int var = 0; var < top; ++var)
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
        delayMili(1);

    }
    QPen red(Qt::red);
    red.setWidth(1);
    forestScene->addLine(-20,0,-1,0,red);
}

/*
 * This functions will add all lines to the scene for relations.
 * This will go through the entire forest and add the lines left
 *  aligned, starting on the 0 axis.
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
    //int startLine = 0;
    int endLine = 0;

    int posYlevel = 0;
    int negYlevel = 0;

    int top = forest->size() - 1;
    int middle = top/2;

    drawRulers(forestScene, penHeight
               , top, 1);

    for (int var = 0; var < middle; ++var)
    {
        //Positive lines
        endLine = (forest->at(middle + var + 1)) / reductionValue;
        //startLine = 0 - endLine;

        posYlevel = ((var*2)+1) *penHeight;

        forestScene->addLine(0//startLine
                             , posYlevel
                             , endLine
                             , posYlevel
                             , green);

        //Negitive lines

        //endLine = ( (forest->at(middle - var -1))/2 ) / reductionValue;
        //startLine = 0 - endLine;

        endLine = forest->at(middle - var -1) / reductionValue;

        negYlevel = ((var*2)+0) * penHeight;

        forestScene->addLine(0 //startLine
                             , negYlevel
                             , endLine
                             , negYlevel
                             , darkGreen);
        delayMili(1);
    }
    QPen red(Qt::red);
    red.setWidth(1);
    forestScene->addLine(-20,0,-1,0,red);

}

void drawRulers(QGraphicsScene *forestScene
                , int penHeight, int forestsize
                , int setOrRelation)
{
    QPen black(Qt::black);
    black.setWidth(1);

    //vertical line
    forestScene->addLine(-20
                         , -20
                         , -20
                         , penHeight*(forestsize -1) //-20
                         , black);
    //horizontal line
    int x = 300;
    forestScene->addLine(-20
                         , -20
                         , 20 + x
                         , -20
                         , black);

    if(setOrRelation == 0)
    {
        //add code for marks on ruler.

        int height = 0;
        //while
    }else
    {


    }
}



/*
 * This functions chooses the default pen height
 *  for use later when adding lines to a scene.
 */
int generatePenHeight(int size)
{
    if(size <20)
    {
        return 15;
    }else if(size >= 20 && size <= 100)
    {
        return 10;
    }else
    {
        return 1;
    }
}

/*
 * This will update the forest scenes.
 *
 */
void drawUpdates(QVector<int> *&forest
                 , QVector<int> &isInForestUpdates
                 , QVector<int> &ForestUpdates
                 , QGraphicsScene *&forestScene
                 , int fPenHeight
                 , int fbase
                 , int fHorizontalReductionValue)
{
    QPen green(Qt::green);
    QPen darkgreen(Qt::darkGreen);

    QGraphicsItem *line = NULL;
    QGraphicsLineItem *lineI = NULL;
    qreal x = 0;
    qreal y = 0;

    //This will process the forest updates.
    //It will decide wheather it is a set or relation
    //then complete the updates to the individual lines.

    //set
    if(fbase > 0)
    {
        while (!ForestUpdates.isEmpty())
        {
            int indexOfForest = ForestUpdates.takeFirst();
            isInForestUpdates.replace(indexOfForest, 0);
            y = indexOfForest * fPenHeight;
            line = forestScene->itemAt(x, y, QTransform());
            lineI = qgraphicsitem_cast<QGraphicsLineItem*>(line);

            int endline = forest->at(indexOfForest);
            endline = endline/fHorizontalReductionValue;

            //If you try to grab a 0 to 0 line then the
            //pointer will be null. This checks to see if
            //you tried to pull an empty line from the
            //scene. If you did then it will create a new
            //line to add to the appropriate spot.
            if(lineI == NULL)
            {
                green.setWidth(fPenHeight);
                forestScene->addLine(0//startline
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
        }
    //Relation
    }else
    {
        while (!ForestUpdates.isEmpty())
        {
            int indexOfForest = ForestUpdates.takeFirst();
            isInForestUpdates.replace(indexOfForest, 0);
            //positive levels
            int middle = forest->size()/2;
            if(indexOfForest > middle)
            {
                int indexOnScene =
                        (indexOfForest - middle) * 2;

                y = indexOnScene * fPenHeight;
                line = forestScene->itemAt(x, y, QTransform());
                forestScene->removeItem(line);
                delete line;
                delayMili(1);

                int endline = forest->at(indexOfForest);
                endline = endline/fHorizontalReductionValue;

                green.setWidth(fPenHeight);
                forestScene->addLine(0//startline
                                      , y
                                      , endline
                                      , y
                                      , green
                                      );

            }else//negitive levels
            {
                int indexOnScene = (forest->size() - 1)
                        - (indexOfForest * 2) -1;
                y = indexOnScene * fPenHeight;
                line = forestScene->itemAt(x, y, QTransform());
                forestScene->removeItem(line);
                delete line;
                delayMili(1);

                int endline = forest->at(indexOfForest);
                endline = endline/fHorizontalReductionValue;

                darkgreen.setWidth(fPenHeight);
                forestScene->addLine(0 //startline
                                      , y
                                      , endline
                                      , y
                                      , darkgreen
                                      );
            }
        }
    }
}

//void drawUpdates(QVector<int> *&forest1
//                 , QVector<int> *&forest2
//                 , QVector<int> &isInForest1Updates
//                 , QVector<int> &isInForest2Updates
//                 , QVector<int> &Forest1Updates
//                 , QVector<int> &Forest2Updates
//                 , QGraphicsScene *&forest1Scene
//                 , QGraphicsScene *&forest2Scene
//                 , int f1PenHeight
//                 , int f2PenHeight
//                 , int f1base
//                 , int f2base
//                 , int f1HorizontalReductionValue
//                 , int f2HorizontalReductionValue)
//{
//    QPen green(Qt::green);
//    QPen darkgreen(Qt::darkGreen);

//    QGraphicsItem *line = NULL;
//    QGraphicsLineItem *lineI = NULL;
//    qreal x = 0;
//    qreal y = 0;

//    //This will process th forest 1 updates.
//    //It will decide wheather it is a set or relation
//    //then complete the updates to the individual lines.

//    //set
//    if(f1base > 0)
//    {
//        while (!Forest1Updates.isEmpty())
//        {
//            int indexOfForest = Forest1Updates.takeFirst();
//            isInForest1Updates.replace(indexOfForest, 0);
//            y = indexOfForest * f1PenHeight;
//            line = forest1Scene->itemAt(x, y, QTransform());
//            lineI = qgraphicsitem_cast<QGraphicsLineItem*>(line);

//            int endline = forest1->at(indexOfForest);
//            endline = endline/f1HorizontalReductionValue;

//            //If you try to grab a 0 to 0 line then the
//            //pointer will be null. This checks to see if
//            //you tried to pull an empty line from the
//            //scene. If you did then it will create a new
//            //line to add to the appropriate spot.
//            if(lineI == NULL)
//            {
//                green.setWidth(f1PenHeight);
//                forest1Scene->addLine(0//startline
//                                      , y
//                                      , endline
//                                      , y
//                                      , green
//                                      );
//            }else//This changes the positionof the
//                //end of the line object.
//            {
//                lineI->setLine(0 //startline
//                               , y
//                               , endline
//                               , y);
//            }
//        }
//    //Relation
//    }else
//    {
//        while (!Forest1Updates.isEmpty())
//        {
//            int indexOfForest = Forest1Updates.takeFirst();
//            isInForest1Updates.replace(indexOfForest, 0);
//            //positive levels
//            int middle = forest1->size()/2;
//            if(indexOfForest > middle)
//            {
//                int indexOnScene = (indexOfForest - middle) * 2;

//                y = indexOnScene * f1PenHeight;
//                line = forest1Scene->itemAt(x, y, QTransform());
//                forest1Scene->removeItem(line);
//                delete line;

//                int endline = forest1->at(indexOfForest);
//                endline = endline/f1HorizontalReductionValue;

//                green.setWidth(f1PenHeight);
//                forest1Scene->addLine(0//startline
//                                      , y
//                                      , endline
//                                      , y
//                                      , green
//                                      );

//            }else//negitive levels
//            {
//                int indexOnScene = (forest1->size() - 1)
//                        - (indexOfForest * 2) -1;
//                y = indexOnScene * f1PenHeight;
//                line = forest1Scene->itemAt(x, y, QTransform());
//                forest1Scene->removeItem(line);
//                delete line;

//                int endline = forest1->at(indexOfForest);
//                endline = endline/f1HorizontalReductionValue;

//                darkgreen.setWidth(f1PenHeight);
//                forest1Scene->addLine(0 //startline
//                                      , y
//                                      , endline
//                                      , y
//                                      , darkgreen
//                                      );
//            }
//        }
//    }

//    //This will process th forest 2 updates.
//    //It will decide wheather it is a set or relation
//    //then complete the updates to the individual lines.

//    //Set
//    if(f2base > 0)
//    {
//        while (!Forest2Updates.isEmpty())
//        {
//            int indexOfForest = Forest2Updates.takeFirst();
//            isInForest2Updates.replace(indexOfForest, 0);
//            y = indexOfForest * f2PenHeight;
//            line = forest2Scene->itemAt(x, y, QTransform());
//            lineI = qgraphicsitem_cast<QGraphicsLineItem*>(line);
//            //forest2Scene->removeItem(line);
//            //delete line;

//            int endline = forest2->at(indexOfForest);
//            endline = endline/f2HorizontalReductionValue;

//            //If you try to grab a 0 to 0 line then the
//            //pointer will be null. This checks to see if
//            //you tried to pull an empty line from the
//            //scene. If you did then it will create a new
//            //line to add to the appropriate spot.
//            if(lineI == NULL)
//            {
//                green.setWidth(f2PenHeight);
//                forest2Scene->addLine(0//startline
//                                      , y
//                                      , endline
//                                      , y
//                                      , green
//                                      );
//            }else//This changes the positionof the
//                //end of the line object.
//            {
//                lineI->setLine(0 //startline
//                               , y
//                               , endline
//                               , y);
//            }


////            green.setWidth(f2PenHeight);
////            forest2Scene->addLine(0 //startline
////                                  , y
////                                  , endline
////                                  , y
////                                  , green
////                                  );
//        }

//    //Relation
//    }else
//    {
//        while (!Forest2Updates.isEmpty())
//        {
//            int indexOfForest = Forest2Updates.takeFirst();
//            isInForest2Updates.replace(indexOfForest, 0);

//            //positive levels
//            int middle = forest2->size()/2;
//            if(indexOfForest > middle)
//            {
//                int indexOnScene = (indexOfForest - middle) * 2;

//                y = indexOnScene * f2PenHeight;
//                line = forest2Scene->itemAt(x, y, QTransform());
//                forest2Scene->removeItem(line);
//                delete line;

//                int endline = forest2->at(indexOfForest);
//                endline = endline/f2HorizontalReductionValue;
//                //int startline = 0 - endline;

//                green.setWidth(f2PenHeight);
//                forest2Scene->addLine(0 //startline
//                                      , y
//                                      , endline
//                                      , y
//                                      , green
//                                      );

//            }else//negitive levels
//            {
//                int indexOnScene = (forest2->size() - 1)
//                        - (indexOfForest * 2) -1;
//                y = indexOnScene * f2PenHeight;
//                line = forest2Scene->itemAt(x, y, QTransform());
//                forest2Scene->removeItem(line);
//                delete line;

//                int endline = forest2->at(indexOfForest);
//                endline = endline/f2HorizontalReductionValue;
//                //int startline = 0 - endline;

//                darkgreen.setWidth(f2PenHeight);
//                forest2Scene->addLine(0 //startline
//                                      , y
//                                      , endline
//                                      , y
//                                      , darkgreen
//                                      );
//            }
//        }
//    }
//}
