#include <QtGui>

#include "mainwindow.h"
#include "log.h"
#include "params.h"
#include "global.h"
//#include "field.h"
#include <QPixmap>
#include <QScrollArea>

#ifdef linux
#include <QTcpServer>
#else
#include <QTcpServer.h>
#endif

LOG_USE();

extern Params *parm;	// I don't believe this is the right way, but it works
QPainter painter;

//QMyGroupBox::QMyGroupBox(QWidget *parent) : QGroupBox(parent)

//class QMyMdiArea : public QScrollArea
//{
//public:
//    QMyMdiArea(QWidget *parent = 0)
//        :
//            QScrollArea(parent)    //, m_pixmap("logo.jpg")

/*
QMyMdiArea::QMyMdiArea(QWidget *parent) : QMdiArea(parent)
    {}
void QMyMdiArea::paintEvent(QPaintEvent *event)
{
        QMdiArea::paintEvent(event);

        QPainter painter(viewport());

        // Calculate the logo position - the bottom right corner of the mdi area.
//        int x = width() - m_pixmap.width();
//        int y = height() - m_pixmap.height();
//        painter.drawPixmap(x, y, m_pixmap);
        QBrush brush;
        brush.setColor(Qt::green);
        painter.setBrush(brush);
        painter.drawEllipse(20,20,10,10);
        painter.end();
        printf("did paint\n");
};
*/

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_comb_DRUG_A_currentIndexChanged(int index)
{
    if (comb_DRUG_A->currentIndex() == -1) return;
    QString filename = comb_DRUG_A->currentText();
    QString drugname = filename;
    QString cell_line;

    if (!cbox_USE_DRUG_A->isChecked()) return;
    extractDrugname(&drugname, &cell_line);
    text_DRUG_A_NAME->setText(drugname);
    radioButton_drugA->setText(drugname);
    if (radioButton_drugA->isChecked()) {
        LOG_QMSG("readDrugParams: " + filename);
        readDrugParams(0, filename);
        LOG_QMSG("populateDrugTable: drug A");
        populateDrugTable(0);
    } else if (radioButton_drugB->isChecked()) {
        LOG_QMSG("populateDrugTable: drug B");
        populateDrugTable(1);
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_comb_DRUG_B_currentIndexChanged(int index)
{
    if (comb_DRUG_B->currentIndex() == -1) return;
    QString filename = comb_DRUG_B->currentText();
    QString drugname = filename;
    QString cell_line;

    if (!cbox_USE_DRUG_B->isChecked()) return;
    extractDrugname(&drugname, &cell_line);
    text_DRUG_B_NAME->setText(drugname);
    radioButton_drugB->setText(drugname);
    if (radioButton_drugB->isChecked()) {
        LOG_QMSG("readDrugParams: " + filename);
        readDrugParams(1, filename);
        LOG_QMSG("populateDrugTable: drug B");
        populateDrugTable(1);
    } else if (radioButton_drugA->isChecked()) {
        LOG_QMSG("populateDrugTable: drug A");
        populateDrugTable(0);
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::extractDrugname(QString *drugname, QString *cell_line)
{
    int i;
    i = drugname->indexOf('.');
    *drugname = drugname->left(i);
    i = drugname->indexOf('_');
    *cell_line = drugname->mid(i+1,99);
    *drugname = drugname->left(i);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_USE_DRUG_A_toggled(bool checked)
{
    QString cell_line;
    LOG_MSG("cbox_use_DRUG_A toggled");
    comb_DRUG_A->setEnabled(checked);
    text_DRUG_A_NAME->setEnabled(checked);
    radioButton_drugA->setEnabled(checked);
    comb_DRUG_A->setCurrentIndex(-1);
    if (checked) {
        radioButton_drugA->setChecked(true);
        /*
        QString drugname = comb_DRUG_A->currentText();
        extractDrugname(&drugname,&cell_line);
        text_DRUG_A_NAME->setText(drugname);
        radioButton_drugA->setText(drugname);
        */
        on_comb_DRUG_A_currentIndexChanged(0);
    } else {
        radioButton_drugA->setChecked(false);
        if (cbox_USE_DRUG_B->isChecked()) {
            radioButton_drugB->click();
            radioButton_drugB->setChecked(true);
        }
        text_DRUG_A_NAME->setText("");
        radioButton_drugA->setText("Drug A");
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_USE_DRUG_B_toggled(bool checked)
{
    QString cell_line;
    LOG_MSG("cbox_use_DRUG_B toggled");
    comb_DRUG_B->setEnabled(checked);
    text_DRUG_B_NAME->setEnabled(checked);
    radioButton_drugB->setEnabled(checked);
    comb_DRUG_B->setCurrentIndex(-1);
    if (checked) {
        radioButton_drugB->setChecked(true);
        /*
        QString drugname = comb_DRUG_B->currentText();
        extractDrugname(&drugname,&cell_line);
        text_DRUG_B_NAME->setText(drugname);
        radioButton_drugB->setText(drugname);
        */
        on_comb_DRUG_B_currentIndexChanged(1);
    } else {
        radioButton_drugB->setChecked(false);
        if (cbox_USE_DRUG_A->isChecked()) {
            LOG_MSG("Check rb_drugA");
            radioButton_drugA->click();
            radioButton_drugA->setChecked(true);
        }
        text_DRUG_B_NAME->setText("");
        radioButton_drugB->setText("Drug B");
    }
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::on_checkBox_CELLDISPLAY_1_toggled(bool display)
{
//    vtk->display_celltype[1] = display;
//    vtk->renderCells();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::on_checkBox_CELLDISPLAY_2_toggled(bool display)
{
//    vtk->display_celltype[2] = display;
//    vtk->renderCells();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::on_comboBox_CELLCOLOUR_1_currentIndexChanged(int index)
{
    QColor qcolor;
    qcolor = comboColour[index];
    Global::celltype_colour[1] = qcolor;    // bad to duplicate this!
    vtk->celltype_colour[1] = qcolor;
    vtk->renderCells();
    sprintf(msg,"changed celltype_colour[1]: index: %d r,g,b: %d %d %d",index,qcolor.red(),qcolor.green(),qcolor.blue());
    LOG_MSG(msg);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::on_comboBox_CELLCOLOUR_2_currentIndexChanged(int index)
{
    QColor qcolor;
    qcolor = comboColour[index];
    Global::celltype_colour[2] = qcolor;    // bad to duplicate this!
    vtk->celltype_colour[2] = qcolor;
    vtk->renderCells();
}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void MainWindow::on_checkBox_FACS_PLOT_toggled(bool checked)
{
    line_FACS_INTERVAL->setEnabled(checked);
    if (!checked) {
        line_FACS_INTERVAL->setText("0");
    }
}

//-------------------------------------------------------------
// Switches to the FACS screen
//-------------------------------------------------------------
void MainWindow::on_action_FACS_triggered()
{
    stackedWidget->setCurrentIndex(4);
    action_outputs->setEnabled(true);
    action_inputs->setEnabled(true);
//    action_VTK->setEnabled(true);
    action_FACS->setEnabled(false);
    Global::showingFACS = true;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_line_CELLPERCENT_1_textEdited(QString pc1_str)
{
    double pc1 = pc1_str.toDouble();
    double pc2 = 100 - pc1;
    QString pc2_str = QString::number(pc2);
    line_CELLPERCENT_2->setText(pc2_str);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_line_CELLPERCENT_2_textEdited(QString pc2_str)
{
    double pc2 = pc2_str.toDouble();
    double pc1 = 100 - pc2;
    QString pc1_str = QString::number(pc1);
    line_CELLPERCENT_1->setText(pc1_str);
}


//------------------------------------------------------------------------------------------------------
// This should be used for any radioButtonGroups for model input parameters
// The thing tp remember about radiobuttons is that ->setChecked(false) programmatically does nothing
// when the radiobutton is one of a groupBox collection.
// This is because in general there may be more than two members of the group, and therefore it is
// impossible to know which one should be set checked.
//------------------------------------------------------------------------------------------------------
void MainWindow::radioButtonChanged(QAbstractButton *b)
{
    QString wtag = b->objectName();
    LOG_QMSG("radioButtonChanged: " + wtag);
    int rbutton_case;
//    if (b->isChecked()) {
        QString ptag = parse_rbutton(wtag,&rbutton_case);
        // Now need to reflect the change in the workingParameterList
        // Need to locate ptag
        LOG_QMSG("ptag: " + ptag);
        wtag = wtag.mid(5);
        for (int k=0; k<nParams; k++) {
            PARAM_SET p = parm->get_param(k);
            if (wtag.compare(p.tag) == 0) {
                LOG_QMSG("found: " + wtag);
                if (ptag.compare("HYPOXIA")==0) {
                    Global::i_hypoxia_cutoff = rbutton_case;
                    QString linetag = "line_HYPOXIA_"+QString::number(rbutton_case);
                    LOG_QMSG("hypoxia tag: " + linetag);
                    QLineEdit *line = findChild<QLineEdit *>(linetag);
                    line_HYPOXIA_THRESHOLD->setText(line->text());
                } else {
                    parm->set_value(k,double(rbutton_case));
                    sprintf(msg,"parm->set_value: %d",rbutton_case);
                    LOG_MSG(msg);
                }
                break;
            }
        }
//        if (radioButton_hypoxia_1->isChecked()) {
//            Global::i_hypoxia_cutoff = 1;
//            line_HYPOXIA_THRESHOLD->setText(line_HYPOXIA_1->text());
//        } else if (radioButton_hypoxia_2->isChecked()) {
//            Global::i_hypoxia_cutoff = 2;
//            line_HYPOXIA_THRESHOLD->setText(line_HYPOXIA_2->text());
//        } else if (radioButton_hypoxia_3->isChecked()) {
//            Global::i_hypoxia_cutoff = 3;
//            line_HYPOXIA_THRESHOLD->setText(line_HYPOXIA_3->text());
//        }

}

void MainWindow::setHypoxiaRadioButton(QString val_str)
{
    if (line_HYPOXIA_1->text().compare(val_str) == 0) {
        rbut_HYPOXIA_1->setChecked(true);
    } else if (line_HYPOXIA_2->text().compare(val_str) == 0) {
        rbut_HYPOXIA_2->setChecked(true);
    } else if (line_HYPOXIA_3->text().compare(val_str) == 0) {
        rbut_HYPOXIA_3->setChecked(true);
    }
}

void MainWindow::buttonClick_cell_constituent(QAbstractButton* button)
{
    LOG_MSG("buttonClick_cell_constituent");
    field->setCellConstituent(button);
}

void MainWindow::buttonClick_field_constituent(QAbstractButton* button)
{
    LOG_MSG("buttonClick_field_constituent");
    field->setFieldConstituent(button);
}

void MainWindow::buttonClick_plane(QAbstractButton* button)
{
    LOG_MSG("buttonClick_plane");
    field->setPlane(button);
}

void MainWindow::buttonClick_canvas(QAbstractButton* button)
{
    LOG_MSG("buttonClick_canvas");
}

void MainWindow::textChanged_fraction(QString text)
{
    LOG_MSG("textChanged_fraction");
    field->setFraction(text);
}

void MainWindow::textEdited_fraction(QString text)
{
    LOG_MSG("textEdited_fraction");
    field->setFraction(text);
}

void MainWindow::onSelectCellConstituent()
{
    if (exthread != NULL) {
        field->selectCellConstituent();
        QString rbname = "rb_cell_constituent_cell" + QString::number(field->cell_constituent);
        QRadioButton *rb = groupBox_cell_constituent->findChild<QRadioButton *>(rbname);
        if (rb) {
            rb->setChecked(true);
        } else {
            LOG_QMSG("onSelectCellConstituent: failed to find rb: " + rbname)
        }
        updateProfilePlots();
    }
}
void MainWindow::onSelectFieldConstituent()
{
    if (exthread != NULL) {
        field->selectFieldConstituent();
        QString rbname = "rb_field_constituent_field" + QString::number(field->field_constituent-1);
        QRadioButton *rb = groupBox_field_constituent->findChild<QRadioButton *>(rbname);
        if (rb) {
            rb->setChecked(true);
        }
        updateProfilePlots();
    }
}

void MainWindow::on_verticalSliderTransparency1_sliderMoved(int position)
{
    vtk->setOpacity(1,position);
}

void MainWindow::on_verticalSliderTransparency2_sliderMoved(int position)
{
    vtk->setOpacity(2,position);
}


//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_checkBox_show_cells_toggled()
{
    int res;
    field->slice_changed = true;
    field->show_cells = checkBox_show_cells->isChecked();
    field->displayField(field->hour,&res);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_checkBox_colony_toggled()
{
    Global::simulate_colony = checkBox_colony->isChecked();
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_lineEdit_colony_days_textChanged()
{
    Global::colony_days = lineEdit_colony_days->text().toDouble();

}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::clickedGraph(QMouseEvent *event)
{
    if (event->button() == Qt::RightButton) {
        LOG_MSG("Right button click");
        QString fileName = QFileDialog::getSaveFileName(this, tr("Select image file"), ".",
            tr("Image files (*.png)"));
        if (fileName.compare("") == 0) {
            return;
        }
        colony_plot->savePng(fileName);
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: on_cbox_SAVE_PROFILE_DATA_toggled(bool checked)
{
    text_SAVE_PROFILE_DATA_FILE_NAME->setEnabled(checked);
    line_SAVE_PROFILE_DATA_INTERVAL->setEnabled(checked);
    line_SAVE_PROFILE_DATA_NUMBER->setEnabled(checked);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: on_cbox_SAVE_SLICE_DATA_toggled(bool checked)
{
    text_SAVE_SLICE_DATA_FILE_NAME->setEnabled(checked);
    line_SAVE_SLICE_DATA_INTERVAL->setEnabled(checked);
    line_SAVE_SLICE_DATA_NUMBER->setEnabled(checked);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow:: on_cbox_SAVE_FACS_DATA_toggled(bool checked)
{
    text_SAVE_FACS_DATA_FILE_NAME->setEnabled(checked);
    line_SAVE_FACS_DATA_INTERVAL->setEnabled(checked);
    line_SAVE_FACS_DATA_NUMBER->setEnabled(checked);
}

//--------------------------------------------------------------------------------------------------------
// Note: with field->displayField, clicking the checkbox can cause a crash.  Not sure why.
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_checkBox_celltypecolours2D_toggled(bool checked)
{
    int res;
    Global::celltypecolours2D = checked;
//    field->slice_changed = true;
//    if (started) field->displayField(hour,&res);
    checkBox_only2colours2D->setEnabled(checked);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_checkBox_only2colours2D_toggled(bool checked)
{
    int res;
    Global::only2colours2D = checked;
//    field->slice_changed = true;
//    if (started) field->displayField(hour,&res);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void MainWindow::on_checkBox_show_bottom_toggled(bool checked)
{
    vtk->show_bottom = checked;
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_pushButton_update_FACS_Histo_clicked()
{
    if (exthread->paused || exthread->stopped) {
        exthread->getFACS();
        emit facs_update();
        emit histo_update();
    }
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_checkBox_volume_scaling_toggled(bool checked)
{
    Global::volume_scaling = checked;
}

//--------------------------------------------------------------------------------------------------------
// For vspheroid this is always true, disabled in the GUI
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_cbox_USE_DIVIDE_TIME_DIST_toggled(bool checked)
{
    groupBox_divisiondistributions->setEnabled(checked);
    line_DVDIVIDE->setEnabled(!checked);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_pushButton_saveHisto_clicked()
{
    saveHistoImage(checkBox_save_histo_data->isChecked());
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::on_pushButton_saveFACS_clicked()
{
    saveFACSImage();
}

//--------------------------------------------------------------------------------------------------------
// This passes the colour scheme selection to the DLL, and to field.
//--------------------------------------------------------------------------------------------------------
void MainWindow::setColourScheme(int icset)
{
    set_colourscheme(icset);
    field->icolourscheme = icset;
    showColours(icset);
}

//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
void MainWindow::showColours(int icset)
{
    int row, col;
    QString valueStr;
    QLabel *label;
    QLineEdit *line;
    COLOUR_SCHEME *cs;

    cs = &field->colScheme[icset];
    for (row = 0; row < 7; row++) {
        label = (QLabel *)(gridLayout_colours->itemAtPosition(row,0)->widget());
        label->setText("");
        for (col=1; col <= 3; col++) {
            line = (QLineEdit *)(gridLayout_colours->itemAtPosition(row,col)->widget());
            line->setText("");
        }
    }
    for (row = 0; row < cs->ncolours; row++) {
        label = (QLabel *)(gridLayout_colours->itemAtPosition(row,0)->widget());
        label->setText(cs->category[row]);
        for (col = 1; col <= 3; col++) {
            line = (QLineEdit *)(gridLayout_colours->itemAtPosition(row,col)->widget());
            switch (col) {
            case 1:
                valueStr = QString::number(cs->color[row].red());
                line->setText(valueStr);
                break;
            case 2:
                valueStr = QString::number(cs->color[row].green());
                line->setText(valueStr);
                break;
            case 3:
                valueStr = QString::number(cs->color[row].blue());
                line->setText(valueStr);
                break;
            }
        }
    }

    if (cs_scene == NULL) {
        cs_layout = new QVBoxLayout();
        cs_scene = new QGraphicsScene();
        cs_view = new QGraphicsView(cs_scene);
        groupBox_show->setLayout(cs_layout);
        cs_layout->addWidget(cs_view);
        connect(lineEdit_sc_0R,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_1R,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_2R,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_3R,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_4R,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_5R,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_6R,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_0G,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_1G,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_2G,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_3G,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_4G,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_5G,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_6G,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_0B,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_1B,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_2B,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_3B,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_4B,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_5B,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
        connect(lineEdit_sc_6B,SIGNAL(textEdited(QString)), this, SLOT(sc_textEdited(QString)));
    }
    cs_scene->clear();
    double dy = 27;
    double w = 25;
    QBrush brush;
    brush.setStyle(Qt::SolidPattern);
    for (row = 0; row < 7; row++) {
        if (row < cs->ncolours) {
            brush.setColor(cs->color[row]);
        } else {
            brush.setColor(QColor(255,255,255));
        }
        cs_scene->addEllipse(5,row*dy,w,w,Qt::NoPen,brush);
    }
    cs_view->show();
}

void MainWindow::sc_textEdited(QString text)
{
    int icset = field->icolourscheme;
    COLOUR_SCHEME *cs = &field->colScheme[icset];
    QObject *obj = sender();
    QString name = obj->objectName();
    QString numstr = name[12];
    int num = numstr.toInt();
    QString colstr = name[13];
    if (colstr.compare("R") == 0) {
        cs->color[num].setRed(text.toInt());
    } else if (colstr.compare("G") == 0) {
        cs->color[num].setGreen(text.toInt());
    } else if (colstr.compare("B") == 0) {
        cs->color[num].setBlue(text.toInt());
    }
    showColours(icset);
}

    /*
    switch (row) {
    case 0:
        label->setText("Growing");
        break;
    case 1:
        label->setText("S-phase");
        break;
    case 2:
        label->setText("M-phase");
        break;
    case 3:
        label->setText("Non-growing");
        break;
    case 4:
        label->setText("Dying (ATP)");
        break;
    case 5:
        label->setText("Dying (treatment)");
        break;
    case 6:
        label->setText("Hypoxic");
        break;
    }
*/
    /*
    brush.setColor(Qt::red);
    cs_scene->addEllipse(5,0,w,w,Qt::NoPen,brush);
    brush.setColor(Qt::green);
    cs_scene->addEllipse(5,dy,w,w,Qt::NoPen,brush);
    brush.setColor(Qt::yellow);
    cs_scene->addEllipse(5,2*dy,w,w,Qt::NoPen,brush);
    brush.setColor(Qt::blue);
    cs_scene->addEllipse(5,3*dy,w,w,Qt::NoPen,brush);
    brush.setColor(Qt::red);
    cs_scene->addEllipse(5,4*dy,w,w,Qt::NoPen,brush);
    brush.setColor(Qt::green);
    cs_scene->addEllipse(5,5*dy,w,w,Qt::NoPen,brush);
    brush.setColor(Qt::yellow);
    cs_scene->addEllipse(5,6*dy,w,w,Qt::NoPen,brush);
    cs_view->show();
    */
