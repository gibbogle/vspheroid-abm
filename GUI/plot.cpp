#include "plot.h"
#include "log.h"

LOG_USE();

#define USE_LEGEND false

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
QString getImageFile()
{
	QString fileName = QFileDialog::getSaveFileName(0,"Select image file", ".", 
		"Image files (*.png *.jpg *.tif *.bmp)");
	return fileName;
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
Plot::Plot(QString aname, QString acasename, QWidget *parent)
	: QwtPlot(parent)
{
	name = aname;
	for (int i=0; i<ncmax; i++) {
		curve[i] = 0;
	}
	ncurves = 0;
	yscale = 0;
    if (acasename.compare("conc")==0) {
        setAxisTitle(QwtPlot::xBottom, "Distance (um)");
    } else if (acasename.compare("vol")==0) {
        setAxisTitle(QwtPlot::xBottom, "Volume fraction");
    } else if (acasename.compare("oxy")==0) {
        setAxisTitle(QwtPlot::xBottom, "O2 level");
    } else {
        setAxisTitle(QwtPlot::xBottom, "Time (hours)");
    }
	if (name.compare("") != 0) {
		curve[0] = new QwtPlotCurve(acasename);
		curve[0]->attach(this);
		ncurves = 1;
		replot();
	}
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
Plot::~Plot()
{
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Plot::mousePressEvent (QMouseEvent *event) {
	event->accept();
	if (event->button() == Qt::RightButton) {
		QString fileName = getImageFile();
		if (fileName.isEmpty()) {
			return;
		}
        QSizeF size(120,120);
        QwtPlotRenderer renderer;
        renderer.renderDocument(this,fileName,size,85);
    }
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Plot::addCurve(QString name)
{
	for (int k=0; k<ncmax; k++) {
		if (curve[k] == 0) {
			curve[k] = new QwtPlotCurve(name);
			curve[k]->attach(this);
			ncurves++;
			replot();
			break;
		}
	}
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Plot::removeCurve(QString name)
{
	for (int i=0; i<ncmax; i++) {
		if (curve[i] != 0) {
			if (name.compare(curve[i]->title().text()) == 0) {
				curve[i]->detach();
				delete curve[i];
				curve[i] = 0;
				ncurves--;
			}
		}
	}
    replot();
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Plot::removeAllCurves()
{
	for (int i=0; i<ncmax; i++) {
		if (curve[i] != 0) {
			curve[i]->detach();
			delete curve[i];
			curve[i] = 0;
			ncurves--;
		}
	}
	replot();
}
	
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Plot::setYScale(double maxval)
{
	yscale = calc_yscale(maxval);
	setAxisScale(QwtPlot::yLeft, 0, yscale, 0);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
double Plot::calc_yscale(double yval)
{
    return 1.3*yval;
//    /*
//    int v = int(1.3*yval);
//    if (v >= 1) {
//        return double(v);
//    } else {
//        return 1.3*yval;
//    }
//    */
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Plot::redraw(double *x, double *y, int n, QString name, QString tag, double fixed_yscale, bool profile)
{
    // Note: Number of pen colors should match ncmax
    QColor pencolor[] = {Qt::black, Qt::red, Qt::blue, Qt::darkGreen, Qt::magenta, Qt::darkCyan };
    QPen *pen = new QPen();
    for (int k=0; k<ncmax; k++) {
        if (curve[k] == 0) continue;
//        if (profile) LOG_QMSG("profile redraw "+name);
        if (name.compare(curve[k]->title().text()) == 0) {
//            LOG_QMSG("redraw " + tag);
            // Just in case someone set ncmax > # of pen colors (currently = 6)
            if (k < 6) {
                pen->setColor(pencolor[k]);
            } else {
                pen->setColor(pencolor[0]);
            }
            curve[k]->setPen(*pen);
            curve[k]->setSamples(x, y, n);
            if (fixed_yscale == 0) {
                double ylast = y[n-1];
                if (ylast > yscale) {
                    yscale = max(yscale,calc_yscale(ylast));
                }
            } else {
                yscale = fixed_yscale;
            }
            setAxisScale(QwtPlot::yLeft, 0, yscale, 0);
            replot();
        }
    }
    delete pen;
}


