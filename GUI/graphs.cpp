#include <qstring.h>
#include "graphs.h"
#include "log.h"

LOG_USE();

Graphs::Graphs()
{
GRAPH_SET tsGraphSet[] = {

    {"nlive",
    "Live Cells",
    "No. of cells",
    "Number of live cells in the blob",
    1, false, 0, 1, 0, TS_TYPE},

    {"nanoxiadead",
    "Anoxia-killed Cells",
    "No. of cells",
     "Total number of cells that have been killed by anoxia",
    2, false, 0, 1, 0, TS_TYPE},

    {"naglucosiadead",
    "Aglucosia-killed Cells",
    "No. of cells",
     "Total number of cells that have been killed by aglucosia",
    3, false, 0, 1, 0, TS_TYPE},

    {"ndrugAdead",
    "DrugA-killed Cells",
    "No. of cells",
     "Total number of cells that have been killed by drugA",
    4, false, 0, 1, 0, TS_TYPE},

    {"ndrugBdead",
    "DrugB-killed Cells",
    "No. of cells",
     "Total number of cells that have been killed by drugB",
    5, false, 0, 1, 0, TS_TYPE},

    {"nradiationdead",
    "Radiation-killed Cells",
    "No. of cells",
     "Total number of cells that have been killed by radiation",
    6, true, 0, 1, 0, TS_TYPE},

    {"nanoxiatagged",
    "Anoxia-tagged Cells",
    "No. of cells",
     "Current number of cells tagged to die by anoxia",
    7, true, 0, 1, 0, TS_TYPE},

    {"naglucosiatagged",
    "Aglucosia-tagged Cells",
    "No. of cells",
     "Current number of cells tagged to die by aglucosia",
    8, true, 0, 1, 0, TS_TYPE},

    {"ndrugAtagged",
    "DrugA-tagged Cells",
    "No. of cells",
     "Current number of cells tagged to die by drugA",
    9, true, 0, 1, 0, TS_TYPE},

    {"ndrugBtagged",
    "DrugB-tagged Cells",
    "No. of cells",
     "Current number of cells tagged to die by drugB",
    10, false, 0, 1, 0, TS_TYPE},

    {"nradiationtagged",
    "Radiation-tagged Cells",
    "No. of cells",
     "Current number of cells tagged to die by radiation",
    11, false, 0, 1, 0, TS_TYPE},

    {"diameter",
    "Spheroid Diameter",
    "Diameter (um)",
     "Spheroid diameter (um)",
    12, true, 0, 1, 0, TS_TYPE},

    {"volume",
    "Spheroid Volume",
    "Volume (mm3)",
     "Spheroid volume (mm3)",
    13, true, 0, 0.001, 0, TS_TYPE},

    {"hypoxicfraction",
    "Hypoxic %",
    "%",
     "Percentage of cells with oxygen level below the specified threshold for hypoxia",
    14, true, 0, 0.1, 0, TS_TYPE},

    {"clonohypoxicfraction",
    "Clonogenic Hypoxic %",
    "%",
     "Percentage of clonogenic cells with oxygen level below the specified threshold for hypoxia",
    15, true, 0, 0.1, 0, TS_TYPE},

    {"growthfraction",
    "Slow-growth Fraction",
    "%",
     "Percentage of cells that are growing at a rate less than the specified fraction of the mean growth rate with no nutrient limits",
    16, false, 0, 0.1, 0, TS_TYPE},

    {"necroticfraction",
    "Necrotic Fraction",
    "%",
     "Percentage of the spheroid that is necrotic = 100*(number of vacant sites)/(number of sites taken up by the spheroid)",
    17, true, 0, 0.1, 0, TS_TYPE},

    {"platingefficiency",
    "Plating Efficiency",
    "%",
     "Plating efficiency = 100*(number of viable cells)/(number of live cells)",
    18, true, 0, 0.1, 0, TS_TYPE},

    {"cellspermm3",
     "Cells/mm3",
     "Density",
     "Number of live cells per mm3 in the blob",
    19, true, 0, 1, 0, TS_TYPE},

    {"mediumoxygen",
    "Ave Medium Oxygen",
    "Concentration",
     "Average concentration of oxygen in the medium (far-field)",
    20, true, 0, 0.001, 0, TS_TYPE},

    {"mediumglucose",
    "Ave Medium Glucose",
    "Concentration",
     "Average concentration of glucose in the medium (far-field)",
    21, true, 0, 0.001, 0, TS_TYPE},

    {"mediumdrugA",
    "Ave Medium Drug A",
    "Concentration",
     "Average concentration of drug A in the medium (far-field)",
    22, true, 0, 0.001, 0, TS_TYPE},

    {"mediumdrugB",
    "Ave Medium Drug B",
    "Concentration",
     "Average concentration of drug B in the medium (far-field)",
    23, true, 0, 0.001, 0, TS_TYPE},

    {"bdryoxygen",
    "Blob Boundary Oxygen",
    "Concentration",
     "Average concentration of oxygen at the blob boundary",
    24, true, 0, 0.001, 0, TS_TYPE},

    {"bdryglucose",
    "Blob Boundary Glucose",
    "Concentration",
     "Average concentration of glucose at the blob boundary",
    25, true, 0, 0.001, 0, TS_TYPE},

    {"bdrydrugA",
    "Blob Boundary Drug A",
    "Concentration",
     "Average concentration of drug A at the blob boundary",
    26, true, 0, 0.001, 0, TS_TYPE},

    {"bdrydrugB",
    "Blob Boundary Drug B",
    "Concentration",
     "Average concentration of drug B at the blob boundary",
    27, true, 0, 0.001, 0, TS_TYPE},

// Extracellular profiles

    {"MULTI",
    "Multi-constituent",
    "",
    "MULTI description",
    MULTI, true, 0, 1, 0, PROF_TYPE},

//    {"CFSE",
//    "CFSE Concentration",
//    "",
//    "CFSE description",
//    CFSE, false, 0, 1, 0, PROF_TYPE},

    {"Oxygen",
    "Oxygen Concentration",
    "",
    "Oxygen description",
    OXYGEN, false, 0, 1, 0, PROF_TYPE},

    {"Glucose",
    "Glucose Concentration",
    "",
    "Glucose description",
    GLUCOSE, false, 0, 1, 0, PROF_TYPE},

    {"Tracer",
    "Tracer Concentration",
    "",
    "Tracer description",
    TRACER, false, 0, 1, 0, PROF_TYPE},

    {"Drug_A",
    "Drug A Concentration",
    "",
    "Drug_A description",
    DRUG_A_PARENT, false, 0, 1, 0, PROF_TYPE},

    {"Drug_A_metab1",
    "Drug A Metabolite 1 Concentration",
    "",
    "Drug_A_metab1 description",
    DRUG_A_METAB_1, false, 0, 1, 0, PROF_TYPE},

    {"Drug_A_metab2",
    "Drug A Metabolite 2 Concentration",
    "",
    "Drug_A_metab2 description",
    DRUG_A_METAB_2, false, 0, 1, 0, PROF_TYPE},

    {"Drug_B",
    "Drug B Concentration",
    "",
    "Drug_B description",
    DRUG_B_PARENT, false, 0, 1, 0, PROF_TYPE},

    {"Drug_B_metab1",
    "Drug B Metabolite 1 Concentration",
    "",
    "Drug_B_metab1 description",
    DRUG_B_METAB_1, false, 0, 1, 0, PROF_TYPE},

    {"Drug_B_metab2",
    "Drug B Metabolite 2 Concentration",
    "",
    "Drug_B_metab2 description",
    DRUG_B_METAB_2, false, 0, 1, 0, PROF_TYPE},

// Intracellular profiles

    {"IC_MULTI",
    "IC Multi-constituent",
    "",
    "IC MULTI description",
    IC_MULTI, true, 0, 1, 0, PROF_TYPE},

    {"IC_Oxygen",
    "IC Oxygen Concentration",
    "",
    "IC Oxygen description",
    IC_OXYGEN, false, 0, 1, 0, PROF_TYPE},

    {"IC_Glucose",
    "IC Glucose Concentration",
    "",
    "IC Glucose description",
    IC_GLUCOSE, false, 0, 1, 0, PROF_TYPE},

    {"IC_Drug_A",
    "IC Drug A Concentration",
    "",
    "Drug_A description",
    IC_DRUG_A_PARENT, false, 0, 1, 0, PROF_TYPE},

    {"IC_Drug_A_metab1",
    "IC Drug A Metabolite 1 Concentration",
    "",
    "Drug_A_metab1 description",
    IC_DRUG_A_METAB_1, false, 0, 1, 0, PROF_TYPE},

    {"IC_Drug_A_metab2",
    "IC Drug A Metabolite 2 Concentration",
    "",
    "Drug_A_metab2 description",
    IC_DRUG_A_METAB_2, false, 0, 1, 0, PROF_TYPE},

    {"IC_Drug_B",
    "IC Drug B Concentration",
    "",
    "Drug_B description",
    IC_DRUG_B_PARENT, false, 0, 1, 0, PROF_TYPE},

    {"Drug_B_metab1",
    "IC Drug B Metabolite 1 Concentration",
    "",
    "Drug_B_metab1 description",
    IC_DRUG_B_METAB_1, false, 0, 1, 0, PROF_TYPE},

    {"Drug_B_metab2",
    "IC Drug B Metabolite 2 Concentration",
    "",
    "Drug_B_metab2 description",
    IC_DRUG_B_METAB_2, false, 0, 1, 0, PROF_TYPE},

    {"IC_CFSE",
    "CFSE Concentration",
    "",
    "IC CFSE description",
    IC_CFSE, false, 0, 1, 0, PROF_TYPE},

    {"IC_growthrate",
    "Growth Rate",
    "",
    "Growth rate description",
    IC_GROWTH_RATE, false, 0, 1, 0, PROF_TYPE},

    {"IC_cellvolume",
    "Cell Volume",
    "",
    "Cell volume description",
    IC_CELL_VOLUME, false, 0, 1, 0, PROF_TYPE},

    {"IC_O2byvolume",
    "Cell O2xVolume",
    "",
    "Cell volume description",
    IC_O2_BY_VOL, false, 0, 1, 0, PROF_TYPE},


// Distributions
//    {"Oxygen_dist",
//    "Oxygen distribution",
//    "Prob",
//    "Oxygen description",
//    OXYGEN, false, 0, 1, 0, DIST_TYPE},

//    {"cellvolume_dist",
//    "Cell volume distribution",
//    "Prob",
//    "Cell volume description",
//    CELL_VOLUME, false, 0, 1, 0, DIST_TYPE}


};

    n_tsGraphs = sizeof(tsGraphSet)/sizeof(GRAPH_SET);
    tsGraphs = new GRAPH_SET[n_tsGraphs];
    for (int i=0; i<n_tsGraphs; i++) {
        tsGraphs[i] = tsGraphSet[i];
    }
    graphList = new GRAPH_SET[maxGraphs];
    nGraphs = maxGraphs;
}


GRAPH_SET Graphs::get_graph(int k)
{
	return graphList[k];
}

int Graphs::get_dataIndex(int k)
{
	return graphList[k].dataIndex;
}

QString Graphs::get_tag(int k)
{
	return graphList[k].tag;
}

QString Graphs::get_title(int k)
{
	return graphList[k].title;
}

QString Graphs::get_yAxisTitle(int k)
{
	return graphList[k].yAxisTitle;
}

QString Graphs::get_description(int k)
{
    return graphList[k].description;
}

double Graphs::get_maxValue(int k) {
	return graphList[k].maxValue;
}

double Graphs::get_scaling(int k) {
	return graphList[k].scaling;
}

double Graphs::get_yscale(int k) {
    return graphList[k].yscale;
}

double Graphs::get_xscale(double xmax) {
    int n = 1;
    for (;;) {
        if (xmax <= n) break;
        n++;
    }
    return double(n);
}

bool Graphs::isActive(int k)
{
	return graphList[k].active;
}

int Graphs::get_type(int k) {
    return graphList[k].type;
}

bool Graphs::isTimeseries(int k)
{
    return (graphList[k].type == TS_TYPE);
}

bool Graphs::isProfile(int k)
{
    return (graphList[k].type == PROF_TYPE);
}

bool Graphs::isDistribution(int k)
{
    return (graphList[k].type == DIST_TYPE);
}

void Graphs::set_maxValue(int k, double v)
{
	graphList[k].maxValue = v;
}

void Graphs::makeGraphList(int non_ts)
{
    int k = maxGraphs;
    int nts = 0;
    for (int i=0; i<n_tsGraphs; i++) {
        if (tsGraphs[i].active) {
            k--;
            graphList[k] = tsGraphs[i];
            nts++;
            if (nts == maxGraphs - non_ts) break;
        }
    }
    int ndummy = maxGraphs - nts - non_ts;
    for (k=0; k<ndummy; k++) {
        graphList[k].tag = "dummy";
        graphList[k].active = false;
        graphList[k].type = TS_TYPE;
        graphList[k].scaling = 1;
    }
    for (k=ndummy; k<ndummy + non_ts; k++) {
        graphList[k].tag = "non_ts";
        graphList[k].active = true;
        graphList[k].type = DIST_TYPE;  //????
        graphList[k].scaling = 1;
    }
    nGraphs = maxGraphs;

    char msg[128];
    sprintf(msg,"nGraphs: %d  non_ts: %d  nts: %d",nGraphs,non_ts,nts);
    LOG_MSG(msg);
//    for (k=0; k<nGraphs; k++) {
//        LOG_QMSG(graphList[k].tag);
//        sprintf(msg,"k: %d scaling: %f",k,graphList[k].scaling);
//        LOG_MSG(msg);
//    }
}

