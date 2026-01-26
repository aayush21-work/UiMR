#include <iostream>
#include <fstream>



int main()
{
  int row, col;
  std::ifstream shap("shape.dat");
  shap >> row >> col ;
  double arr[row][col];

  std::ifstream file("bounds1.dat");
  for(int i=0;i<row;i++)
  {
    for(int j=0;j<col;j++)
    {
      file >> arr[i][j];
    }
  }

  double tempx;
  double tempy;

  for(int i=0;i<row-1;i++)
  {
    tempx=(arr[i][1]+arr[i+1][0])/2.0;
    arr[i][1]=tempx;
    arr[i+1][0]=tempx;

  }


  for(int j=0;j<col-1;j++)
  {
    tempy=(arr[j][2]+arr[j][3])/2.0;
    arr[j][2]=tempy;
    arr[j][3]=tempy;
  }


  std::ofstream out_file("out_file.dat");

  for(int i=0;i<row;i++)
  {
    for(int j=0;j<col;j++)
    {
      out_file << arr[i][j] << " " ;
    }

    out_file << std::endl;
  }



  return 0;

} 
