// Read location of atoms and calculate charge densities on a grid 
#include<iostream>
#include<fstream>
#include<string>
#include <iomanip>
//#include<conio.h>
#include<math.h>
#include<stdio.h>
using namespace std;
//const int Nx=16,Ny=16,Nz=256;
#define I(i,j,k) (k)*Nx*Ny + (j)*Nx +(i) 	//Indexing array

int main(int argc, char **argv)
{   
 if (argc != 3){
    std::cout << argv[0] << " [charge file +] [charge file -] \n";
    return 0;
  }

  int Nf , Nx , Ny , Nz , Na , Npc , Nion , Ncnt , fa;   	

//inputs--------------------------------------------------	
	std::ifstream ifile1;
        ifile1.open("inputs_density.txt");
        ifile1 >> Nf >> Nx >> Ny >> Nz >> Na >> Npc >> Nion >> Ncnt >> fa;
        ifile1.close();
//read files-----------------------------------------------


     std::ofstream ofs("density.txt");   					//output file- choose the name and location of the output file here
     std::ifstream ifs1("coord.crd");   					// file pointer for reading coordinates- select name and location of input file
     std::ifstream ifs1n("coord.crd"); 					//using seperate file pointer for finding box vectors (need this info in advance)- select name and location of the same input file
     std::ifstream ifs2("PC.txt");  						//charge files
     std::ifstream ifs3("TEA.txt");
     std::ifstream ifs4("BF4.txt");
     std::ifstream ifs5(argv[1]);
     std::ifstream ifs6(argv[2]);
//store charges in arrays-------------------------------------------

float q_PC[14]; //13 atoms PC
for(int i=1; i<=13;i++)
ifs2 >> q_PC[i];
//cout<<"\nq[13]="<<q_PC[13]<<"\n";

float q_TEA[30]; //29 atoms in TEA 
for(int j=1; j<=29;j++)
ifs3 >> q_TEA[j];

float q_BF4[6]; //5 atoms  
for(int k=1; k<=5;k++)
	{
ifs4>>q_BF4[k];}

float q_CNTP[361]; //360 atoms defined as double, the chrge values are too long 
for(int l=1; l<=360;l++)
{
ifs5 >> q_CNTP[l];
}

float q_CNTN[361]; //360 atoms defined as double
for(int m=1; m<=360;m++)
{ifs6 >> q_CNTN[m];
}

//Read first(title) line from crdbox-----------------------------------------------------------------------------------
string line;
getline (ifs1,line);//this file pointer to read coordinates
cout << line << '\n';
//Read first(title) line from crdbox
getline (ifs1n,line);//this file pointer to find box vectors
//---------------------------------------------------------------------------------------------------------------------
 
 float X_max,Y_max,Z_max,junk;
 float ax,ay,az,p1,q1,r1,Ax,Ay,Az,x1,y1,z1;
 float x=0,y=0,z=0;

 float *q = new float[Nx*Ny*Nz] ();

//Counter
 int counter=Nf;
   
for(int f=1;f<=Nf;f++) //loop over frames-------------------------------------------------------------------------------
        {
         // volume of each small unit is 1.333 * 1.333 * approx(1.333)
         //loop to find the max of frame
         for(int temp=0; temp<((Na*3));temp++)
                  ifs1n>>junk; //throwing these values
         ifs1n>>X_max;//saving X_max of frame
         ifs1n>>Y_max;//saving Y_max of frame        
         ifs1n>>Z_max;//saving Z_max of frame
         
	 if(f==1) // calculate these only for f=1. Density is calculated at grid points excluding the boundaries at (x=Lx/2), (Y=Ly/2) and (z=Lz/2Y=Ly/2) and (z=Lz/2). Values are copied into these boundaries later (PBC)
         {
	 ax=X_max/(fa*(Nx-1));//factor of 2 because of 1/4 of box
         ay=Y_max/(fa*(Ny-1));
         Ax=1/ax;
         Ay=1/ay;
         x1=(X_max/fa)-(ax)*0.5;
         y1=(Y_max/fa)-(ay)*0.5;
         }
         
	az=Z_max/(Nz-1);       
         Az=1/az;
         z1=(Z_max-(az*0.5));//no factor of 2 here         
         
                  
         for(int pc=1;pc<=Npc;pc++)
                 {for(int at=1;at<=13;at++)
                          {
                                    ifs1>>x;    ifs1>>y;    ifs1>>z;
                                    
                                    if(x>X_max/fa)  // copying all coordinates into 1/4 box
                                    x=x-X_max/fa;
                                    if(y>Y_max/fa)
				    y=y-Y_max/fa;
                          
                                    p1=(x>x1)?(x-x1):x;
                                    q1=(y>y1)?(y-y1):y;
                                    r1=(z>z1)?(z-z1):z;
                                    q[I(int((p1+(ax*0.5))*Ax),int((q1+(ay*0.5))*Ay),int((r1+(az*0.5))*Az))]+=q_PC[at]*Ax*Ay*Az;

                                              
                          }
                 }
                                     
         for(int tea=1;tea<=Nion;tea++)
                 {for(int at=1;at<=29;at++)
                          {
                                          
                                    ifs1>>x;    ifs1>>y;    ifs1>>z;
                                    
			            if(x>X_max/fa)  // copying all coordinates into 1/4 box
                                    x=x-X_max/fa;
                                    if(y>Y_max/fa)
				    y=y-Y_max/fa;
                                    
                                    p1=(x>x1)?(x-x1):x; //
                                    q1=(y>y1)?(y-y1):y;
                                    r1=(z>z1)?(z-z1):z;
                                    q[I(int((p1+(ax*0.5))*Ax),int((q1+(ay*0.5))*Ay),int((r1+(az*0.5))*Az))]+=q_TEA[at]*Ax*Ay*Az;
                          }
                 }  
         for(int bf4=1;bf4<=Nion;bf4++)
                 {for(int at=1;at<=5;at++)
                          {
                                    ifs1>>x;    ifs1>>y;    ifs1>>z;
                                   
				    if(x>X_max/fa)  // copying all coordinates into 1/4 box
                                    x=x-X_max/fa;
                                    if(y>Y_max/fa)
				    y=y-Y_max/fa;
                                    
                                    p1=(x>x1)?(x-x1):x;
                                    q1=(y>y1)?(y-y1):y;
                                    r1=(z>z1)?(z-z1):z;
                                    q[I(int((p1+(ax*0.5))*Ax),int((q1+(ay*0.5))*Ay),int((r1+(az*0.5))*Az))]+=q_BF4[at]*Ax*Ay*Az;
                          }
                 }               
         for(int cntp=1;cntp<=Ncnt;cntp++)
                 {for(int at=1;at<=360;at++)
                          {
                                    ifs1>>x;    ifs1>>y;    ifs1>>z;
                                    
				    if(x>X_max/fa)  // copying all coordinates into 1/4 box
                                    x=x-X_max/fa;
                                    if(y>Y_max/fa)
				    y=y-Y_max/fa;
	
                                    p1=(x>x1)?(x-x1):x;
                                    q1=(y>y1)?(y-y1):y;
                                    r1=(z>z1)?(z-z1):z;
                                    q[I(int((p1+(ax*0.5))*Ax),int((q1+(ay*0.5))*Ay),int((r1+(az*0.5))*Az))]+=q_CNTP[at]*Ax*Ay*Az;
                          }
                 } 
         for(int cntn=1;cntn<=Ncnt;cntn++)
                 {for(int at=1;at<=360;at++)
                          {
                                    ifs1>>x;    ifs1>>y;    ifs1>>z;
                                    if(x>X_max/fa)  // copying all coordinates into 1/4 box
                                    x=x-X_max/fa;
                                    if(y>Y_max/fa)
				    y=y-Y_max/fa;
                                    
                                    p1=(x>x1)?(x-x1):x;
                                    q1=(y>y1)?(y-y1):y;
                                    r1=(z>z1)?(z-z1):z;
                                    q[I(int((p1+(ax*0.5))*Ax),int((q1+(ay*0.5))*Ay),int((r1+(az*0.5))*Az))]+=q_CNTN[at]*Ax*Ay*Az;
                          }
                 }                                 
                   ifs1>>junk;ifs1>>junk;ifs1>>junk;
}

cout<<"\n";
//---------------------------------------------------------------------------------------------------------------------------------------------------------------
//PBC //to copy values into the boundaries
//end points
           q[I(0,0,Nz-1)]=q[I(0,0,0)];
           q[I(0,Ny-1,0)]=q[I(0,0,0)];
           q[I(Nx-1,0,0)]=q[I(0,0,0)];

//edges    
//long       
 int c1,d1,h1,e1,f1;
           for(c1=0;c1<Nz;c1++)
           {
            q[I(0,Ny-1,c1)]=q[I(0,0,c1)];
            q[I(Nx-1,0,c1)]=q[I(0,0,c1)];
            }
//short
            for(d1=0;d1<Nx;d1++)
           {
            q[I(d1,Ny-1,0)]=q[I(d1,0,0)];
            q[I(d1,0,Nz-1)]=q[I(d1,0,0)];
            }
            for(h1=0;h1<Ny;h1++)
           {
            q[I(0,h1,Nz-1)]=q[I(0,h1,0)];
              q[I(Nx-1,h1,0)]=q[I(0,h1,0)];
            }
            
//1st long face
            for(e1=0;e1<Nz;e1++)
            {
                    for(f1=0;f1<Ny;f1++)
                    {         q[I(Nx-1,f1,e1)]=q[I(0,f1,e1)]; 
                    } 
            }         
//2nd long face
            for(c1=0;c1<Nz;c1++)
            {
                    for(d1=0;d1<Nx;d1++)
                    {     
                        q[I(d1,Ny-1,c1)]=q[I(d1,0,c1)];
                    }
            }    
//short face   
               for(c1=0;c1<Ny;c1++)
               {
                     for(d1=0;d1<Nx;d1++)
                     {
                      q[I(d1,c1,Nz-1)]=q[I(d1,c1,0)];
                      }
               }                                        
//-----------------------------------------------------------------------------------------------------------------------------------pbc conditions end
   for(int k=0;k<Nz;k++)
   {
           for(int j=0;j<Ny;j++)   
           {
                   for(int i=0;i<Nx;i++)
                   {
                 
                   ofs<<q[I(i,j,k)]<<"  ";                  
                   }
                   ofs<<"\n";
                 
   
        }       
  }  
ifs1n.close();
ifs1.close();
ifs2.close();
ifs3.close();
ifs4.close();
ifs5.close();
ifs6.close();
ofs.close();

        
//getch();
delete[] q;
return 0;
}
