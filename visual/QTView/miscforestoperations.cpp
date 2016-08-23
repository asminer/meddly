#include "miscforestoperations.h"
#include "drawfunctions.h"
#include "timefunctions.h"

#include <QDebug>
#include <QMessageBox>


QVector<int> *setupRelation(QStringList &list)
{
    QVector<int> *forest = new QVector<int>();

    int value = 0;
    QString temp = NULL;

    for (int var = 5; var < list.size(); ++var)
    {
        temp = list.at(var);
        value = temp.toInt();
        forest->append(value);
    }
    return forest;
}

QVector<int> *setupSet(QStringList &list)
{
    QVector<int> *forest = new QVector<int>();

    int value = 0;
    QString temp = NULL;

    //forest->append(1);
    for (int var = 5; var < list.size(); ++var)
    {
        temp = list.at(var);
        value = temp.toInt();
        forest->append(value);
    }
    return forest;
}

QVector<int> *setupForest1(Parser *parser
                           , int &f1base
                           , int &f1PenHeight
                           , QString &f1Name)
{
    //This is to skip the comment line.
    parser->readLine();

    //The first F line
    QString line = parser->readLine();

    //Split the F line
    // "[,\[ ]"
    QRegExp rx("[,\[ \\]]");// match a comma or a space or square brackets
    QStringList list = line.split(rx, QString::SkipEmptyParts);

    //Setting the forest number and name
    f1Name = "Forest 1 " + list.at(2);

    QString n = list.at(3);
    f1base = n.toInt();
    if(f1base > 0)
    {
        //Call generatePenHeight
        QString d = list.at(4);
        f1PenHeight = generatePenHeight(d.toInt());

        return setupSet(list);
    }else
    {
        //Call generatePenHeight
        QString d = list.at(4);
        f1PenHeight = generatePenHeight(d.toInt() * 2);

        return setupRelation(list);
    }
}

QVector<int> *setupForest2(Parser *parser
                            , int &f2base
                            , int &f2PenHeight
                            , QString &f2Name)
{
    //The second F line
    QString line = parser->readLine();

    //Split the F line
    // "[,\[ ]"
    QRegExp rx("[,\[ \\]]");// match a comma or a space or square brackets
    QStringList list = line.split(rx, QString::SkipEmptyParts);

    //Setting the forest number and name
    f2Name = "Forest 2 " + list.at(2);
    QString n = list.at(3);
    f2base = n.toInt();
    if(f2base > 0)
    {
        //Call generatePenHeight
        QString d = list.at(4);
        f2PenHeight = generatePenHeight(d.toInt());

        return setupSet(list);
    }else
    {
        //Call generatePenHeight
        QString d = list.at(4);
        f2PenHeight = generatePenHeight(d.toInt() * 2);

        return setupRelation(list);
    }
}

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
                 , int &oldMicroseconds
                 , QString f1Name
                 , QString f2Name
                 , bool &endOfFile
                 , int f1HorizontalReductionValue
                 , int f2HorizontalReductionValue
                 , int speedFactorValue)
{

    bool continueLoop = true;
    while(continueLoop)
    {
        QString line = parser->readLine();

        if (line.length() <1)
        {
            //QMessageBox::information(NULL,"File Name", "End of file");
            endOfFile = true;

            //Empty the last of the updates
            drawUpdates(forest1
                        , forest2
                        , isInForest1Updates
                        , isInForest2Updates
                        , Forest1Updates
                        , Forest2Updates
                        , forest1Scene
                        , forest2Scene
                        , f1PenHeight
                        , f2PenHeight
                        , f1base
                        , f2base
                        , f1HorizontalReductionValue
                        , f2HorizontalReductionValue);
            delete parser;
            delete forest1;
            delete forest2;

            break;
        }
        QChar t = line.at(0);

        int u = t.toLatin1();

        switch (u)
        {
        //This is the p line
        case 112:
        {
            QChar fNumber = line.at(2);
            if(fNumber == '1')
            {
                QString phase = line.remove(0, 3);
                labelForest1->setText(f1Name+ phase);
            }

            if(fNumber == '2')
            {
                QString phase = line.remove(0, 3);
                labelForest2->setText(f2Name + phase);
            }
            continueLoop = false;

            break;
        }
        //This is the 'a' line
        case 97:
        {
            QStringList updates = line.split(" ", QString::SkipEmptyParts);
            QString whichForest = NULL;
            QString sexpectedLevel = NULL;
            QString schangeAmount = NULL;

            int indexFromLog = 0;
            for (int var = 1; var < updates.size(); ++var)
            {
                whichForest = updates.at(var);

                ++var;
                sexpectedLevel = updates.at(var);
                indexFromLog = sexpectedLevel.toInt();

                ++var;
                schangeAmount = updates.at(var);
                int changeAmount = schangeAmount.toInt();

                if(whichForest == "1")
                {
                    int indexOfForest = indexFromLog - f1base;
                    if(isInForest1Updates.at(indexOfForest) == 0)
                    {
                        Forest1Updates.append(indexOfForest);
                        isInForest1Updates.replace(indexOfForest, 1);
                    }

                    changeAmount = changeAmount + forest1->at(indexOfForest);
                    forest1->replace(indexOfForest, changeAmount);
                }else
                {
                    int indexOfForest = indexFromLog - f2base;
                    if(isInForest2Updates.at(indexOfForest) == 0)
                    {
                        Forest2Updates.append(indexOfForest);
                        isInForest2Updates.replace(indexOfForest, 1);
                    }

                    changeAmount = changeAmount + forest2->at(indexOfForest);
                    forest2->replace(indexOfForest, changeAmount);
                }

            }
            //continueLoop = false;
            break;
        }

            //This is the t line
        case 116:
        {
            QStringList times = line.split(" ", QString::SkipEmptyParts);

            QString sseconds = times.at(1);
            int newSeconds = sseconds.toInt();
            QString sMicro = times.at(2);
            int newMicroseconds = sMicro.toInt();

            int tempSeconds = newSeconds - oldSeconds;
            int tempMicroseconds = newMicroseconds - oldMicroseconds;

            if(tempSeconds > 0 || tempMicroseconds >= (10000 / speedFactorValue))
            {
                oldSeconds = newSeconds;

                delaySec(tempSeconds);

                //Calculate speed factor slowdown
                //tempMiliseconds = tempMiliseconds*speedFactorValue;
                oldMicroseconds = newMicroseconds;

                //convert mili to micro
                int tempMiliseconds = tempMicroseconds/1000;

                delayMili(tempMiliseconds*speedFactorValue);

                drawUpdates(forest1
                            , forest2
                            , isInForest1Updates
                            , isInForest2Updates
                            , Forest1Updates
                            , Forest2Updates
                            , forest1Scene
                            , forest2Scene
                            , f1PenHeight
                            , f2PenHeight
                            , f1base
                            , f2base
                            , f1HorizontalReductionValue
                            , f2HorizontalReductionValue);

            }
            continueLoop = false;
            break;
        }
        default:
            QMessageBox::information(NULL,"File Name", line);
            continue;
        }//end switch
    }//end while loop
}

void redrawForest(QGraphicsScene *forestScene
                  , QVector<int> *forest, int base
                  , int penHeight, int reductionValue)
{
    forestScene->clear();

    if (base > 0)
    {
        drawSet(forestScene, forest
                , penHeight, reductionValue);
    }else
    {
        drawRelation(forestScene, forest
                     , penHeight, reductionValue);
    }

}
