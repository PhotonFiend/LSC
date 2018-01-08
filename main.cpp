//All measurements to be taken in um.
//Last modified 08/01/2017 by Christopher Richards
//Email : cjr945@uowmail.edu.au

#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


using namespace std;



typedef vector <double> record_t;
typedef vector <record_t> data_t;

//Update!
//Window is now a user defined option.
int UpperWindow =700;
int LowerWindow =650;

int geochoice;
double SpotRadius, DyeConc, AbsWeight1, AbsWeight2, StepSize = 10, ScalingFactor, FluxFactor;
double MaximumAbsWavelength, PowerIn = 0, PowerOut = 0, QuantYield = 0, QuantYieldCheck = 0;
double AbsPeakVal = 0, ExtCoeff,MaxAbsCoeff;
int solarphotoncount = 0;
double X, Y, Z, R;
const time_t ctt = time(0);
bool activephoton;
string ResultsTableName;

data_t SolarEmsData;
data_t DyeAbsData;
data_t DyeEmsData;

//read csv and split into fields (first is wavelength)

istream& operator >> ( istream& ins, record_t& record )
{
    // make sure that the returned record contains only the stuff we read now
    record.clear();

    // read the entire line into a string (a CSV record is terminated by a newline)
    string line;
    getline( ins, line );

    // now we'll use a stringstream to separate the fields out of the line
    stringstream ss( line );
    string field;
    while (getline( ss, field, ',' ))
    {
        // for each field we wish to convert it to a double
        // (since we require that the CSV contains nothing but floating-point values)
        stringstream fs( field );
        double f = 0.0;  // (default value is 0.0)
        fs >> f;

        // add the newly-converted field to the end of the record
        record.push_back( f );
    }

    // Now we have read a single line, converted into a list of fields, converted the fields
    // from strings to doubles, and stored the results in the argument record, so
    // we just return the argument stream as required for this kind of input overload function.
    return ins;
}

//-----------------------------------------------------------------------------
// Let's likewise overload the stream input operator to read a list of CSV records.

istream& operator >> ( istream& ins, data_t& data )
{
    // make sure that the returned data only contains the CSV data we read here
    data.clear();

    // For every record we can read from the file, append it to our resulting data
    record_t record;
    while (ins >> record)
    {
        data.push_back( record );
    }

    // Again, return the argument stream as required for this kind of input stream overload.
    return ins;
}


void WriteFinalPos(double x, double y, double z, double photonWavelength)
{
    ofstream PosFile;
    PosFile.open(ResultsTableName.c_str(), ios::app);
    PosFile << x << "," << y << "," << z << "," << photonWavelength << endl;
    PosFile.close();
}


void WriteAbsValue(double z, double AbsWavelength)
{
    ofstream AbsFile;
    AbsFile.open ("AbsEvent.txt", ios::app);
    AbsFile << z << "," << AbsWavelength << endl;
    AbsFile.close();

}

void WriteSolInput(double SolarWavelength)
{
    ofstream SolFile;
    SolFile.open ("SolarInput.txt", ios::app);
    SolFile <<SolarWavelength << endl;
    SolFile.close();

}

void WriteEmsValue(double x, double y, double z, double AbsWavelength, double EmsWavelength)
{
    ofstream EmsFile;
    EmsFile.open ("EmsEvent.txt", ios::app);
    EmsFile <<x <<"," << y << "," << z << "," << AbsWavelength<< ","<< EmsWavelength << endl;
    EmsFile.close();

}

void userinput()
{
    //simple user input sub for geometery conditions

    bool validchoice = false;

    while (validchoice == false)
    {
        cout << "What is the main body geometry?:\n1. Box\n2. Cylinder" <<endl;
        cin >> geochoice;
        cout << endl;
        if (geochoice == 1) //Box
        {
            double x,y,z;
            cout << "The following dimensions are in mm.\n";
            cout << "Enter X:";
            cin >> x;
            cout << endl;
            cout << "Enter Y:";
            cin >> y;
            cout << endl;
            cout << "Enter Z:";
            cin >> z;
            cout << endl;
            X=x*1000; //Converts units from mm to um
            Y=y*1000;
            Z=z*1000;
            validchoice = true;
        }

        else if (geochoice == 2) //cylinder
        {
            double r,z;
            cout << "The following dimensions are in mm.\n";
            cout << "Enter Radius(R):";
            cin >> r;
            cout << endl;
            cout << "Enter Depth(Z):";
            cin >> z;
            cout << endl;
            validchoice = true;
            R=r*1000;
            Z=z*1000;
        }
        else
        {
            validchoice = false;
            cout << "Not a valid choice, please try again." << endl;
        }
    }

    cout << "\nExtinction Coefficient(/M/cm): ";
    cin >> ExtCoeff;

    cout << "\nDye Concentration(uM): ";
    cin >> DyeConc;

    cout << "\nPlease enter the scaling factor(1 to 27): ";
    cin >> ScalingFactor;
    if(ScalingFactor > 5) cout << "/nWarning: Scaling Factors over 5 may take a long time." <<endl;

    cout <<"\nUpper window on wavelength limit: ";
    cin >> UpperWindow;

    cout << "\nLower window on wavelength limit: ";
    cin >> LowerWindow;

    cout << "\nQuantum Yield (as percentage %): ";
    cin >> QuantYield;

}

double TotalSolarPower()
{
    double ScaledFlux, UserArea, SolarPowerTotal =0;
    int n = 0, j =0 ;
    //for calc of full solar spectrum (outside of what is used in experiment)

    UserArea = X*Y/1000000000000; //area in m^2

    FluxFactor = UserArea * pow(10,ScalingFactor-43);

    for (unsigned n = 0; n < SolarEmsData.size(); n++)
    {
        ScaledFlux = SolarEmsData[n][1] * FluxFactor;
        for( j=0; j< ScaledFlux; j++)
        {
            SolarPowerTotal += 1240/(SolarEmsData[n][0]);
        }
    }
    cout << SolarPowerTotal << "eV" <<endl;
    return (SolarPowerTotal);
}

void loadSolarEms()
{
    // Here is the file containing the data. Read it into memory.
    ifstream infile( "solarems.csv" );
    infile >> SolarEmsData;

    // Show message if something went wrong.
    if (!infile.eof())
    {
        cout << "Error reading solar spectrum\n";

    }

    infile.close();

    cout << "Successfully read " << SolarEmsData.size() << " wavelengths.\n";

    unsigned max_record_size = 0;
    for (unsigned n = 0; n < SolarEmsData.size(); n++)
        if (max_record_size < SolarEmsData[ n ].size())
            max_record_size = SolarEmsData[ n ].size();
    cout << "The largest Ems record has " << max_record_size << " fields.\n";

    return ;
}

void loadDyeAbs()
{
    // Here is the file containing the data. Read it into memory.
    ifstream infile( "dyeabs.csv" );
    infile >> DyeAbsData;

    // Show message if something went wrong.
    if (!infile.eof())
    {
        cout << "Error reading dye absorption spectrum\n";

    }

    infile.close();

    cout << "Successfully read " << DyeAbsData.size() << " wavelengths.\n";

    unsigned max_record_size = 0;
    for (unsigned n = 0; n < DyeAbsData.size(); n++)
        if (max_record_size < DyeAbsData[ n ].size())
            max_record_size = DyeAbsData[ n ].size();
    cout << "The largest Abs record has " << max_record_size << " fields.\n";

    return ;
}

void loadDyeEms()
{
    // Here is the file containing the data. Read it into data.
    ifstream infile( "dyeems.csv" );
    infile >> DyeEmsData;

    // Complain if something went wrong.
    if (!infile.eof())
    {
        cout << "Error reading dye emission spectrum\n";
    }

    infile.close();

    cout << "Successfully read " << DyeAbsData.size() << " wavelengths.\n";

    unsigned max_record_size = 0;
    for (unsigned n = 0; n < DyeAbsData.size(); n++)
        if (max_record_size < DyeAbsData[ n ].size())
            max_record_size = DyeAbsData[ n ].size();
    cout << "The largest Ems record has " << max_record_size << " fields.\n";

    return;
}

void FindMaxAbs()
{
    int n = 0;
    for (unsigned n = 0; n < DyeAbsData.size(); n++)
    {
        if(DyeAbsData[n][1]>=MaxAbsCoeff) MaxAbsCoeff = DyeAbsData[n][1];
    }

}

double GenDyeEms()
{
    //Generate an emitted photon wavelength from the weighted emission
    //spectrum of a dye
    double RandomValue, GenPhoton;
    int n = 0,fraction, datasetsize;
    bool ValueFound;

    ValueFound = false;
    datasetsize = DyeEmsData.size();
    RandomValue = (double(rand())/double(RAND_MAX))*(DyeEmsData[datasetsize - 1][2]);

    fraction = 2;
    while (ValueFound != true)
    {
        if(n<datasetsize-1)
        {
            if(DyeEmsData[n][2]>RandomValue) ValueFound = true;
            else n++;
        }
        else ValueFound = true;
    }
    //GenPhoton = DyeEmsData[n][0];
    //cout << DyeEmsData[n][0] << endl;
    /*
        while (ValueFound != true)
        {
            if (n <= (datasetsize-2))
            {
                if(2*fraction>datasetsize) fraction = datasetsize/4;

                if ((RandomValue <= DyeEmsData[n+1][2]) & (RandomValue > DyeEmsData[n][2]))
                {
                    GenPhoton = DyeEmsData[n][0];
                    ValueFound = true;
                }
                else
                {
                    if (DyeEmsData[n][2] >= RandomValue)
                    {
                        n -= datasetsize/(2*fraction);
                        if (n<=0) n=datasetsize/2;
                        fraction += 1;
                    }
                    else if (DyeEmsData[n][2] <= RandomValue)
                    {
                        n += datasetsize/(2*fraction);
                        fraction += 1;
                            if(n>datasetsize-1) //fix search if exceeds array size
                            {
                                n = datasetsize-2;
                            }
                    }
                }
            }
            else
            {
                n = datasetsize-2;
            }

        }
    */

    return DyeEmsData[n][0] ;
}

void GiveAbsWeight(double wavelength)
{
    //compares the wavelength to the dyeabs dataset for
    //the weighted chance for the photon to be absorbed and
    //re emitted. (Binary search)
    // eg find the weighted value for A for a given wavelength
    // store the values n and n+1 to use as a range.
    // check if the random value is within this range.
    //Add scaling for concentration. Linear Response
    int n,fraction, datasetsize;
    double StdDyeConc, cellvalue, AbsMultiplier; //uM
    bool ValueFound;

    StdDyeConc = 1000000;
    ValueFound = false;
    datasetsize = DyeAbsData.size();
    n = (datasetsize)/2;
    fraction = 2;

    if(DyeConc <= 0)
    {
        AbsWeight1 = 0;
        AbsWeight2 = 0.0000001;
    }

    else
    {


        while (ValueFound != true)
        {

            if ((wavelength <= DyeAbsData[n+1][0]) & (wavelength > DyeAbsData[n][0]))
            {
                //finds the weighted range for the dye absorption.
                //also applies linear transformation of range
                //according to the user specified dye concentration.

                AbsWeight1 = DyeAbsData[n][2];
                AbsWeight2 = DyeAbsData[n+1][2];
                AbsMultiplier = (1-pow(10,(-1*(DyeConc/StdDyeConc)*ExtCoeff*0.001)))*100;

                //AbsWeight2 = (1-AbsMultiplier)*(pow(AbsWeight2-AbsWeight1,2))/(AbsMultiplier*DyeAbsData[datasetsize-1][2]);
                //AbsWeight1 = 0;

                AbsWeight2 = AbsWeight1 + (AbsWeight2- AbsWeight1)* AbsMultiplier;

                ValueFound = true;

            }
            else
            {

                cellvalue=DyeAbsData[n][0];
                //cout << cellvalue << "nm" << endl;
                if (DyeAbsData[n][0] >= wavelength)
                {
                    n -= datasetsize/(2*fraction);
                    fraction += 1;
                    if(n<0) n=0;
                }
                else if (DyeAbsData[n][0] <= wavelength)
                {
                    n += datasetsize/(2*fraction);
                    fraction += 1;
                    if(n>datasetsize) //fix search if exceeds array size
                    {
                        n = datasetsize-2;
                    }
                }
            }

        }

    }
}

bool insidedetector(double x, double y, double z)
{
    //checks if photon position is now inside detector
    // the detector in this case is the edges of the geometry
    //two cases 1. box 2. cylinder

    bool insidebox = false;
    if(geochoice == 1) //in case where prism, use these conditions
    {
        if (((x<-X/2) || (x>X/2)) || ((y<-Y/2) || (y>Y/2)))
        {
            if ((z < 0) || (z > Z)) insidebox = false; //redundancy
            else insidebox = true;
        }
        else
        {
            insidebox = false;
        }
    }
    if(geochoice == 2) //in case where cylinder, use these conditions
    {
        if(sqrt(pow(x,2)+pow(y,2))>R) //circle about origin with radius R
        {
            if((z < 0) || (z > Z)) insidebox = false;//redundancy
            else insidebox = true;
        }
    }

    return insidebox;

}

bool insidebody(double x, double y, double z)
{
    //For regular prisms only the Z value determines the boundary
    //as all edge faces are coupled
    bool insidebox;

    if ((z>Z) || (z<0))
    {
        insidebox = false;
    }
    else insidebox = true;

    return insidebox;
}


void Prism()
{
    int i,j,detectorcount = 0, photoncount= 0, photoncountWindow = 0, ScaledFlux, TotalNumSteps = 0;
    double theta,upsi,x,y,z,reflect,abschance, UserArea, emitchance;
    double photonwavelength, dotProduct, PhotonVector[3] = {0,0,0}, PlaneNorm[3] = {0,0,0};
    double detPhoton, detPlane, AngleCritical;
    bool detected, atboundary, photonReflected;
    detectorcount = 0;
    i=0, j=0;
    int quenches = 0;

    //Standardise the AM1.5 solar spectrum over our given area from user input
    UserArea = X*Y/1000000000000; //area in m^2

    FluxFactor = UserArea * pow(10,ScalingFactor-43);

    photoncount=0;
    for(i=0; i<SolarEmsData.size(); i++)
    {
        //Number of photons per wavelength scaled by user and geometric conditions
        ScaledFlux = SolarEmsData[i][1] * FluxFactor;
        for (j=0; j<ScaledFlux; j++)
        {
            photoncount++;

            //generate starting position based on geometry

            x=(double(rand())/double(RAND_MAX))*X - X/2;//make X and Y value
            y=(double(rand())/double(RAND_MAX))*Y - Y/2;


            //starting face at z=0, moving in +z direction
            z=0;
            activephoton=true;
            //generate photon and the abs of photon in dye
            photonwavelength = SolarEmsData[i][0];

            //add photon energy to power calculations
            PowerIn += 1240/(photonwavelength);
//WriteSolInput(photonwavelength);
            //accounts for solar photons within relevant band gaps
            if((photonwavelength >LowerWindow) & (photonwavelength <UpperWindow)) photoncountWindow++;

            //check to see if solar photon is in AbsData wavelengths
            if(photonwavelength < MaximumAbsWavelength)
            {
                GiveAbsWeight(photonwavelength);

                //first direction vector will not be straight back out of the box
                upsi = 0;
                theta = 90;
                //reset some other variables
                photonReflected = false;
                PlaneNorm[0] = 0;
                PlaneNorm[1] = 0;
                PlaneNorm[2] = 0;

                while(activephoton != false)
                {
                    TotalNumSteps++;
                    //random variables for abs
                    abschance = (double(rand())/RAND_MAX)*(DyeAbsData[DyeAbsData.size()-1][2]);

                    //photon steps by 10um in direction
                    //step method depends on if photon has reflected
                    // off surface.
                                        if (photonReflected == true)
                    {
                        x +=  PhotonVector[0];
                        y += PhotonVector[1];
                        z += PhotonVector[2];

                    }
                    else
                    {
                        x += StepSize*cos(theta)*sin(upsi);
                        y += StepSize*sin(theta)*sin(upsi);
                        z += StepSize*cos(upsi);
                    }

                    detected = insidedetector(x,y,z);
                    atboundary = insidebody(x,y,z);

                    if (detected == true) //checks if photon is in UC volume
                    {
                        activephoton = false; //flag set false
                        WriteFinalPos(x,y,z,photonwavelength);//save photon info
                        PowerOut += 1240/(photonwavelength);
                        detectorcount++;
                    }
                    else if (abschance > AbsWeight1 & abschance < AbsWeight2)
                    {
                        //photon absorbed by dye and remitted at new wavelength in
                        //dye emission spectrum. Also scattered.
                        //is there a chance for non-radiative relaxation over certain wavelength?

                        //WriteAbsValue(z,photonwavelength);

                        QuantYieldCheck = 100*(double(rand())/double(RAND_MAX));
                        if(QuantYieldCheck>QuantYield)
                        {
                            activephoton = false;
                            quenches++;
                        }
                        else
                        {
                            double prevWavelength = photonwavelength;

                            photonwavelength = GenDyeEms();
                            if(photonwavelength > DyeEmsData[DyeEmsData.size()-2][0])
                            {
                                activephoton = false;
                            }
                            //option for only emitting photons above current wavelength
                            /*while(photonwavelength < prevWavelength)
                            {
                                photonwavelength = GenDyeEms();
                            }*/

                            //WriteEmsValue(x,y,z,prevWavelength,photonwavelength);



                            GiveAbsWeight(photonwavelength);

                            upsi = (double(rand())/RAND_MAX)*2*3.1415;
                            theta = acos((double(rand())/RAND_MAX)*2-1);
                        }
                        //WriteEmsValue(photonwavelength);
                    }
                    else if (atboundary == false) //photon has reached the boundary of body
                    {
                        //check to see if photon at coupled edge
                        //it should be picked up by insideDetector() before this
                        //lets be safe anyway.
                        if (((x<-X/2) || (x>X/2)) || ((y<-Y/2) || (y>Y/2)))
                        {
                            activephoton = false;
                            WriteFinalPos(x,y,z,photonwavelength);//save photon info
                            PowerOut += (300000000)*(6.63*pow(10,-34))/(photonwavelength*pow(10,-9));
                            detectorcount++;
                        }
                        else
                        {
                            if (z<0) PlaneNorm[2] = 1;
                            else if (z>Z) PlaneNorm[2] = -1;
                            else if (x>X/2) PlaneNorm[0] = -1; //these other reflections are redundant.
                            else if (x<-X/2) PlaneNorm[0] = 1; //Perhaps useful for non-edge coupled systems.
                            else if (y>Y/2) PlaneNorm[1] = -1; //also would need to clear array after each loop.
                            else if (y<-Y/2) PlaneNorm[1] = 1;

                            if (photonReflected == false) //if the photon hasn't been reflected generate array values
                            {

                            PhotonVector[0] = -StepSize*cos(theta)*sin(upsi);
                            PhotonVector[1] = -StepSize*sin(theta)*sin(upsi);
                            PhotonVector[2] = -StepSize*cos(upsi);
                            }
                            dotProduct = PhotonVector[0]*PlaneNorm[0] + PhotonVector[1]*PlaneNorm[1] + PhotonVector[2]*PlaneNorm[2];
                            detPhoton = sqrt(pow(PhotonVector[0],2)+pow(PhotonVector[1],2)+pow(PhotonVector[2],2)); //ARGHH square roots!
                            detPlane = sqrt(pow(PlaneNorm[0],2)+pow(PlaneNorm[1],2)+pow(PlaneNorm[2],2)); //why do you do this to me???
                            AngleCritical = 90-180*acos((dotProduct)/(detPhoton*detPlane))/3.1415;

                            if(AngleCritical < 43.0 && AngleCritical >-43.0)
                            {
                                //steps photon back inside the body before assigning new direction vectors
                                //prevents out of boundary error

                            x -=  PhotonVector[0];
                            y -=  PhotonVector[1];
                            z -=  PhotonVector[2];


                                photonReflected = true;
                                //calculate angle between plane and photon vector and reflect by equal degrees

                                //PhotonVector[0] = -2*(dotProduct)*PlaneNorm[0] - PhotonVector[0];
                                //PhotonVector[1] = -2*(dotProduct)*PlaneNorm[1] - PhotonVector[1];
                                if (z<0) z=0;
                                if (z>Z) z=Z;
                                PhotonVector[2] = -PhotonVector[2];



                            }

                            else activephoton = false;

                        }


                    }

                }
            }

        }

        cout << "\r" << 100*i/SolarEmsData.size() << "%";
    }

    cout << "\nThe total number of photons at the edges was: " << detectorcount << " fleebs" << endl;
    //output results and parameters to file
    ofstream myfile;
    myfile.open ("Results.txt", ios::app);
    myfile << asctime(localtime(&ctt)) << endl;
    myfile << "Test conditions\n===============\nX\t\tY\t\tZ\tConc.\tStpSize" << endl;
    myfile << X << "\t"<< Y << "\t"<< Z << "\t" << DyeConc << "\t" << StepSize << endl;
    myfile << "\nNumber of photons: " << photoncount << endl;
    myfile << "\nPhoton Count at edges: " << detectorcount << endl;
    myfile << "\nPhotons/cm^2 top surface: " << (photoncount/(X*Y))*100000000 << endl;
    myfile << "\nPhotons/cm^2 edges: " << (detectorcount/(2*((Z*Y)+(Z*X))))*100000000 << endl;
    myfile << "\nPhoton Flux at edges: " << (detectorcount/(2*((Z*Y)+(Z*X))))/(photoncount/(X*Y))*100 << "% of AM1.5 incident intensity" << endl;
    myfile << "\nPhoton Flux at edges: " << (detectorcount/(2*((Z*Y)+(Z*X))))/(photoncountWindow/(X*Y))*100 << "% of AM1.5 incident intensity between" << endl;
    myfile << "\nPower input top: " << PowerIn << endl;
    myfile << "\nPower output edges: " << PowerOut << endl;
    myfile << "\nPhoton input power/area(inside simulation window): " << (PowerIn/(X*Y))*100000000 << endl;
    myfile << "\nPhoton input power/area(total solar): " << (TotalSolarPower()/(X*Y))*100000000 << endl;
    myfile << "\nPhoton output power/area: " << (PowerOut/(2*((Z*Y)+(Z*X))))*100000000 <<endl;
    myfile << "\nPower Ratio Side/Top (Geometric): "<<100*(((PowerOut/(2*((Z*Y)+(Z*X))))*100000000)/((TotalSolarPower()/(X*Y))*100000000)) << "%"<<endl;
    myfile << "\nPower Ratio Side/Top: " <<100*PowerOut/PowerIn << "%" <<endl;
    myfile << "\nAverage Photon Path Length: " << double(TotalNumSteps*StepSize/photoncount) << " um" <<endl;
    myfile << "===========================\n\n" <<endl;

    myfile.close();
    cout << quenches <<endl;
    cin.ignore();

}

void Cylinder()
{
    int i,j,detectorcount = 0, photoncount= 0, photoncountWindow = 0, ScaledFlux, TotalNumSteps = 0;
    double theta,upsi,x,y,z,reflect,abschance, UserArea, emitchance;
    double photonwavelength;
    bool activephoton, detected = false, atboundary;
    detectorcount = 0;
    i=0, j=0;
    //Standardise the AM1.5 solar spectrum over our given area from user input

    UserArea = R*R*3.1415/1000000000000;
    FluxFactor = UserArea * pow(10,ScalingFactor);

    photoncount=0;
    for(i=0; i<SolarEmsData.size(); i++)
    {
        //Number of photons per wavelength scaled by user and geometric conditions
        ScaledFlux = SolarEmsData[i][1] * FluxFactor;
        for (j=0; j<ScaledFlux; j++)
        {
            photoncount++;

            //generate starting position based on geometry
            x=(double(rand())/double(RAND_MAX))*R*2 - R;//make X and Y value
            y=(double(rand())/double(RAND_MAX))*R*2 - R;
            while(sqrt(x*x+y*y)>R) //make sure X and Y fall within radius
            {
                x=(double(rand())/double(RAND_MAX))*R*2 - R; //regenerate until within radius
                y=(double(rand())/double(RAND_MAX))*R*2 - R;
            }

            //starting face at z=0, moving in +z direction
            z=0;
            activephoton=true;
            //generate photon and the abs of photon in dye
            photonwavelength = SolarEmsData[i][0];

            //accounts for solar photons within relevant band gaps
            if((photonwavelength >LowerWindow) & (photonwavelength <UpperWindow)) photoncountWindow++;

            //check to see if solar photon is in AbsData wavelengths
            if(photonwavelength < MaximumAbsWavelength)
            {
                GiveAbsWeight(photonwavelength);

                //first direction vector will not be straight back out of the box
                upsi = 0;
                theta = 90;

                while(activephoton != false)
                {
                    TotalNumSteps++;
                    //random variables for abs
                    abschance = (double(rand())/RAND_MAX)*(DyeAbsData[DyeAbsData.size()-1][2]);
                    reflect = (double(rand())/RAND_MAX);
                    //photon steps by 1um in direction
                    x += StepSize*cos(theta)*sin(upsi);
                    y += StepSize*sin(theta)*sin(upsi);
                    z += StepSize*cos(upsi);

                    detected = insidedetector(x,y,z);
                    atboundary = insidebody(x,y,z);

                    if (detected == true) //checks if photon is in UC volume
                    {
                        activephoton = false; //flag set false
                        WriteFinalPos(x,y,z,photonwavelength);//save photon info
                        detectorcount++;
                    }
                    else if (abschance > AbsWeight1 & abschance < AbsWeight2)
                    {
                        //photon absorbed by dye and remitted at new wavelength in
                        //dye emission spectrum. Also scattered.
                        //is there a chance for non-radiative relaxation over certain wavelength?

                        //WriteAbsValue(photonwavelength);
                        QuantYieldCheck = 100*(double(rand())/double(RAND_MAX));
                        if(QuantYieldCheck>QuantYield)
                        {
                            activephoton = false;
                        }
                        else
                        {

                            photonwavelength = GenDyeEms();
                            GiveAbsWeight(photonwavelength);

                            upsi = (double(rand())/RAND_MAX)*2*3.1415;
                            theta = acos((double(rand())/RAND_MAX)*2-1);
                        }
                        //WriteEmsValue(photonwavelength);
                    }
                    else if (atboundary == false) //photon has reached the boundary of body



                        if (reflect < 0.75) //chance photon is reflected internally
                        {
                            if(z>Z)
                            {
                                activephoton = false;
                            }
                            else
                            {
                                upsi = (double(rand())/RAND_MAX)*2*3.1415;
                                theta = acos((double(rand())/RAND_MAX)*2-1);
                            }
                        }
                        else
                        {
                            //photon escaped. byebye photon.
                            activephoton = false;
                        }
                }

            }
        }

    }

    cout << "\r" << 100*i/SolarEmsData.size() << "%";


cout << "\nThe total number of photons at the edges was: " << detectorcount << " fleebs" << endl;
//output results and parameters to file
ofstream myfile;
myfile.open ("Results.txt", ios::app);
myfile << asctime(localtime(&ctt)) << endl;
myfile << "Test conditions\n===============\nR\tZ\tSpotSize\tConc.\tStpSize" << endl;
myfile << R << "\t"<< Z << "\t" <<SpotRadius << "\t\t" << DyeConc << "\t" << StepSize << endl;
myfile << "\nNumber of photons: " << photoncount << endl;
myfile << "\nPhoton Count at edges: " << detectorcount << endl;
myfile << "\nPhotons/cm^2 top surface: " << (photoncount/(R*R*3.1415))*100000000 << endl;
myfile << "\nPhotons/cm^2 edges: " << (detectorcount/((2*3.1415*R*Z)))*100000000 << endl;
myfile << "\nPhoton Flux at edges: " << (detectorcount/((2*3.1415*R*Z)))/(photoncount/(R*R*3.1415))*100 << "% of AM1.5 incident intensity" << endl;
//the following line is soft-coded with the solar photon count (650-700nm) for debugging purposes.
//only for scaling factor = 2
myfile << "\nPhoton Flux at edges: " << (detectorcount/(2*(R*R*3.1415)))/(photoncountWindow/(R*R*3.1415))*100 << "% of AM1.5 incident intensity between" << endl;
myfile << "\nAverage Photon Path Length: " << double(TotalNumSteps*StepSize/photoncount) << " um" <<endl;
myfile << "===========================\n\n" <<endl;

myfile.close();
cin.ignore();

}

void outputDataFileSetup()
{
    //Generate output filename based on user input and variables
    cout << "\nOutput Results table: ";
    stringstream TempString;
    TempString << "ResultTable_";
    TempString << DyeConc;
    TempString << "uM_";
    TempString <<  X/1000;
    TempString << "x";
    TempString << Y/1000;
    TempString << "x";
    TempString << Z/1000;
    TempString << "mm.txt";
    ResultsTableName += TempString.str();
    cout << ResultsTableName <<endl;
    cout << endl;

}

int main()
{
    srand(time(0));
    cout << "      MonteCarlo UC\n";
    cout << "-------------------------\n";
    cout << "Calculating Paths...\nApproximating Universe...\nClearing Sentience Cache..." <<endl;
    //load the required data into memory
    //Future versions would allow the user to choose csv's manually
    loadSolarEms();
    loadDyeAbs();
    loadDyeEms();
    FindMaxAbs();



    MaximumAbsWavelength = DyeAbsData[DyeAbsData.size()-1][0]; // The max wavelength from given data

    //let the user pick the boundary conditions
    userinput();
    outputDataFileSetup();
    if (geochoice == 1) Prism(); //area in m^2
    if (geochoice == 2) Cylinder();


    return 0;

}


