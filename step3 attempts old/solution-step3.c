// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2019.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2019 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>


double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;


/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-4) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  
  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}



/**
 * This is the only operation you are allowed to change in the assignment.
 */
 //used to check is the previous iteration had a collision
void updateBody() {
  const int bucketNumber = 10;
  int** buckets = new int*[bucketNumber](); //initialise a pointer to pointers for all [bucketNumber] buckets
  for(int i=0; i<bucketNumber; i++){
    buckets[i] = new int[NumberOfBodies](); //each bucket can hold up to [NumberOfBodies] values
  }
  int* bucketPointers = new int[bucketNumber]();
  
  if (!(t>0)){
    for(int i=0; i<NumberOfBodies;i++){
      buckets[0][bucketPointers[0]] = i;
      bucketPointers[0] += 1;
    }
  }else{
    double vBucket = maxV/bucketNumber;
    for (int i=0; i<NumberOfBodies; i++){
      double velocity = std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] );
      if(velocity <= vBucket){
        buckets[0][bucketPointers[0]] = i;
        bucketPointers[0] += 1;
      } else if (velocity <= vBucket * 2){
        buckets[1][bucketPointers[1]] = i;
        bucketPointers[1] += 1;
      } else if (velocity <= vBucket * 3){
        buckets[2][bucketPointers[2]] = i;
        bucketPointers[2] += 1;
      } else if (velocity <= vBucket * 4){
        buckets[3][bucketPointers[3]] = i;
        bucketPointers[3] += 1;
      } else if (velocity <= vBucket * 5){
        buckets[4][bucketPointers[4]] = i;
        bucketPointers[4] += 1;
      } else if (velocity <= vBucket * 6){
        buckets[5][bucketPointers[5]] = i;
        bucketPointers[5] += 1;
      } else if (velocity <= vBucket * 7){
        buckets[6][bucketPointers[6]] = i;
        bucketPointers[6] += 1;
      } else if (velocity <= vBucket * 8){
        buckets[7][bucketPointers[7]] = i;
        bucketPointers[7] += 1;
      } else if (velocity <= vBucket * 9){
        buckets[8][bucketPointers[8]] = i;
        bucketPointers[8] += 1;
      } else {
        buckets[9][bucketPointers[9]] = i;
        bucketPointers[9] += 1;
      }
    }
  }
  
  static bool prevCol = false; //used to check is the previous iteration had a collision
  const double threshold = 0.01*0.01;
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();
  // force0 = force along x direction
  // force1 = force along y direction
  // force2 = force along z direction
  double* force0 = new double[NumberOfBodies]();
  double* force1 = new double[NumberOfBodies]();
  double* force2 = new double[NumberOfBodies]();
  //store the previous forces for adams-bashforth
  static double* prevf0 = new double[NumberOfBodies]();
  static double* prevf1 = new double[NumberOfBodies]();
  static double* prevf2 = new double[NumberOfBodies]();
  

  for(int k=0; k<bucketNumber; k++){
    int timesteps = pow(2,k);
    double alteredTime = timeStepSize/timesteps;
    for(int step=0; step<timesteps; step++){
      for (int jj=0; jj<bucketPointers[k]; jj++){ //loop through all particles in bucket
        int j = buckets[k][jj];
        if(j != -1){
          force0[j] = 0.0;
          force1[j] = 0.0;
          force2[j] = 0.0;
          for (int i=0; i<NumberOfBodies; i++) { //for each other particle do
            if (i != j){
              //storing calculations that are use multiple times to avoid redundancy
              const double dist0 = x[i][0]-x[j][0], dist1 = x[i][1]-x[j][1], dist2 = x[i][2]-x[j][2];
              const double distance = sqrt(dist0*dist0 + dist1*dist1 + dist2*dist2);

              // x,y,z forces acting on particles
              
              const double a = mass[i]*mass[j] /(distance * distance * distance) ; 
              //effect of other paritcle on current particle
              const double f0 = dist0 * a, f1 = dist1 * a, f2 = dist2 * a;
              force0[j] += f0 ;
              force1[j] += f1 ;
              force2[j] += f2 ;
              minDx = std::min( minDx,distance );
            }
          }
          x[j][0] = x[j][0] + alteredTime * v[j][0];
          x[j][1] = x[j][1] + alteredTime * v[j][1];
          x[j][2] = x[j][2] + alteredTime * v[j][2];
          //Using Adams-Bashforth for velocity to increase accuracy
          //not used if the first iteration or if a collision happened in the previous iteration (as the prevForce value will be wrong)
          if(t<0 && !prevCol){

            v[j][0] = v[j][0] + alteredTime*(1.5 * force0[j]/mass[j] - 0.5 * prevf0[j]/mass[j]);
            v[j][1] = v[j][1] + alteredTime*(1.5 * force1[j]/mass[j] - 0.5 * prevf1[j]/mass[j]);
            v[j][2] = v[j][2] + alteredTime*(1.5 * force2[j]/mass[j] - 0.5 * prevf2[j]/mass[j]);

          }else{

            v[j][0] = v[j][0] + alteredTime * force0[j] / mass[j];
            v[j][1] = v[j][1] + alteredTime * force1[j] / mass[j];
            v[j][2] = v[j][2] + alteredTime * force2[j] / mass[j];
            prevCol = false;
          }
          maxV = std::max(maxV,std::sqrt( v[j][0]*v[j][0] + v[j][1]*v[j][1] + v[j][2]*v[j][2] ));
          prevf0[j] = force0[j];
          prevf1[j] = force1[j];
          prevf2[j] = force2[j];
        }
      }
      //collision checker
      for (int j=0; j<NumberOfBodies; j++){
        for(int i=j+1; i<NumberOfBodies; i++){
          const double dist0 = x[i][0]-x[j][0], dist1 = x[i][1]-x[j][1], dist2 = x[i][2]-x[j][2];
          const double squareDist = dist0*dist0 + dist1*dist1 + dist2*dist2;
          if (squareDist <= threshold){ 
            const double newMass = mass[i]+mass[j];
            for(int n=0; n<3; n++){ 
              v[j][n] = (mass[i]*v[i][n] + mass[j]*v[j][n])*(1/newMass);
              x[j][n] = (x[i][n] + x[j][n]) * 0.5;
            }
            mass[j] = newMass;

            if(i!=NumberOfBodies-1){//if the collided particle is not the end one, store the end particles data in place of it (as we are no longer using it)
              mass[i] = mass[NumberOfBodies-1];
              x[i] = x[NumberOfBodies-1];
              v[i] = v[NumberOfBodies-1];
            }
            prevCol = true;
            for(int bucket=0; bucket<bucketNumber; bucket++){
              for(int particle=0; particle<bucketPointers[bucket]; particle++){
                if(buckets[bucket][particle] = i){
                  buckets[bucket][particle] = -1;
                } else if (buckets[bucket][particle] = NumberOfBodies -1){
                  buckets[bucket][particle] = i;
                }
              }
            }
            NumberOfBodies -= 1;

          }
        }
      }
    }
  }

  
  //collision checker, using the square of the distance as a threshold (less error is involved when not square rooting)
  
  
  if (NumberOfBodies==1){
    std::cout <<  x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
    t=tFinal+1;   
  }

  t += timeStepSize;

  delete[] force0;
  delete[] force1;
  delete[] force2;
  delete[] buckets;
  delete[] bucketPointers;
}


/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << "  dt              time step size (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t Five body setup" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;

      tPlot += tPlotDelta;
    }
  }

  closeParaviewVideoFile();

  return 0;
}
