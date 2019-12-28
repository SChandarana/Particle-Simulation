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
double smallestMin = std::numeric_limits<double>::max();
double prevmaxV;
int whichcol = 0;
void updateBody() {
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();
  bool collision = false;
  double threshold = 0.01 * 0.01;
  // force0 = force along x direction
  // force1 = force along y direction
  // force2 = force along z direction
  auto* force0 = new double[NumberOfBodies]();
  auto* force1 = new double[NumberOfBodies]();
  auto* force2 = new double[NumberOfBodies]();
  static auto* prevf0 = new double[NumberOfBodies]();
  static auto* prevf1 = new double[NumberOfBodies]();
  static auto* prevf2 = new double[NumberOfBodies]();
  static auto* prevv0 = new double[NumberOfBodies]();
  static auto* prevv1 = new double[NumberOfBodies]();
  static auto* prevv2 = new double[NumberOfBodies]();
  
  for (int j=0; j<NumberOfBodies; j++){ //loop through all particles
    for (int i=j+1; i<NumberOfBodies; i++) { //for each other particle do (calculations for previous particles are already done thus start at j+1)
      const double dist0 = x[i][0]-x[j][0], dist1 = x[i][1]-x[j][1], dist2 = x[i][2]-x[j][2];
      //square rooting introduces more error so storing the squared version for use where possible
      const double squareDist = dist0*dist0 + dist1*dist1 + dist2*dist2;
      const double distance = sqrt(squareDist);
      // if (squareDist <= threshold){
      //   collision = true;
      //   whichcol += 1;
      //   std::cout << "COLLISION" << std::endl;
      //   smallestMin = std::numeric_limits<double>::max();
      // }
      
      // x,y,z forces acting on particles
      // reducing number of redundant calculations by storing this
      const double a = mass[i]*mass[j] /(squareDist * distance) ; 
      //effect of other paritcle on current particle
      const double f0 = dist0 * a, f1 = dist1 * a, f2 = dist2 * a;
      force0[j] += f0 ;
      force1[j] += f1 ;
      force2[j] += f2 ;
      //effect of current particle on other particle
      force0[i] -= f0 ;
      force1[i] -= f1 ;
      force2[i] -= f2 ;
      minDx = std::min( minDx,distance );
      // if (minDx <= 0.0125 && whichcol != 0){
      //   std::cout << minDx << std::endl;
      //   std::cout << prevmaxV << std::endl;
      //   getchar();
      // }
      
    }
  }

  for(int i=0; i<NumberOfBodies; i++){
    
    
    //using Adam-Bashford as outlined in the lectures
    if (t < 0){
      x[i][0] = x[i][0] + timeStepSize * (1.5 * v[i][0] - 0.5 * prevv0[i]);
      x[i][1] = x[i][1] + timeStepSize * (1.5 * v[i][1] - 0.5 * prevv1[i]);
      x[i][2] = x[i][2] + timeStepSize * (1.5 * v[i][2] - 0.5 * prevv2[i]);
      v[i][0] = v[i][0] + timeStepSize*(1.5 * force0[i]/mass[i] - 0.5 * prevf0[i]/mass[i]);
      v[i][1] = v[i][1] + timeStepSize*(1.5 * force1[i]/mass[i] - 0.5 * prevf1[i]/mass[i]);
      v[i][2] = v[i][2] + timeStepSize*(1.5 * force2[i]/mass[i] - 0.5 * prevf2[i]/mass[i]);
    }else{
      x[i][0] = x[i][0] + timeStepSize * v[i][0];
      x[i][1] = x[i][1] + timeStepSize * v[i][1];
      x[i][2] = x[i][2] + timeStepSize * v[i][2];
      v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
      v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
      v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];
    }

    maxV = std::max(maxV,std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] ));
    prevf0[i] = force0[i];
    prevf1[i] = force1[i];
    prevf2[i] = force2[i];
    prevv0[i] = v[i][0];
    prevv1[i] = v[i][1];
    prevv2[i] = v[i][2];
    prevmaxV = maxV;
  }

  //checking if a collision has occurred
  //start at end of list as removing last element when collision occurs
  //  if(collision){
    for(int i=NumberOfBodies-1; i>0; i--){
      for(int j=0; j<i; j++){
        const double dist0 = x[j][0]-x[i][0], dist1 = x[j][1]-x[i][1], dist2 = x[j][2]-x[i][2];
        //square rooting introduces more error so storing the squared version for use where possible
        const double squareDist = dist0*dist0 + dist1*dist1 + dist2*dist2;
        
        if (squareDist <= threshold){
          std::cout << "HERE" << std::endl;
          const double massNew = mass[j]+mass[i];
          const double invMass = 1/massNew;
          const double massJinv = invMass * mass[j], massIinv = invMass * mass[i];
          //storing new particle information in earlier occuring particle(j)
          //changing position of new particle to average of collided
          x[j][0] = (x[i][0] + x[j][0])*0.5; 
          x[j][1] = (x[i][1] + x[j][1])*0.5;
          x[j][2] = (x[i][2] + x[j][2])*0.5;
          //change velocities of new particle
          v[j][1] = massJinv*v[j][1] + massIinv*v[i][1];
          v[j][2] = massJinv*v[j][2] + massIinv*v[i][2];
          v[j][0] = massJinv*v[j][0] + massIinv*v[i][0];
          maxV = std::max(maxV,std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] ));
          
          //replace the collided particle with the end particle and lower the number of particles(effectively deleting the collided particle)
          //by doing this we can account for a >2-way collision
          if(i != NumberOfBodies -1){
            x[i] = x[NumberOfBodies-1];
            v[i] = v[NumberOfBodies-1];
            mass[i] = mass[NumberOfBodies-1];
            
          }
          NumberOfBodies -= 1;
        }  
      }
  }
  //  }
  if (NumberOfBodies == 1){
    //stop if only one particle left
    std::cout << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
    t = tFinal;
  }

  t += timeStepSize;

  delete[] force0;
  delete[] force1;
  delete[] force2;
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
