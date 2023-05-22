#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<sstream>

using namespace std;

#define PI 3.1421592
#define MaxRows 500
#define MaxCols 500

struct grayImage{

    grayImage()
    {
        Rows = Cols = 0;
        Loaded = false;
        for(int r = 0; r< MaxRows; r++)
            for(int c = 0; c< MaxCols; c++)
                Image[r][c] = 0;
    }

    unsigned short setPixel(unsigned short value, int r, int c)
    {
        if( r >= Rows || c >= Cols || r < 0 || c < 0)
            return -1;
        Image[r][c] = value;
        return value;
    }

    int getPixel(int r, int c)
    {
        if( r >= Rows || c >= Cols || r < 0 || c < 0)
            return -1;
        return Image[r][c];
    }

    int setRows(int rows)
    {
        if(rows < 1 || rows > MaxRows)
            return -1;
        Rows = rows;
        return Rows;
    }

    int getRows()
    {
        return Rows;
    }

    int setCols(int cols)
    {
        if(cols < 1 || cols > MaxCols)
            return -1;
        Cols = cols;
        return Cols;
    }

    int getCols()
    {
        return Cols;
    }

    int Save(string File_Name)
    {
	if (!Loaded)
	{
	    return -1;
	}
	if (File_Name.substr(File_Name.size()-4,4)!=".pgm")
	{
	    return 1;
	}
	ofstream FWRITE(File_Name.c_str());
	if (!FWRITE)
	{
    	    return 2;
	}
	FWRITE<<"P2"<<endl<<"# File Created by The Fast and Curious"<<endl<<Cols<<" "<<Rows<<endl<<Maximum<<endl;
	for (int i=0; i<Rows; i++)
	{
	    for (int j=0; j<Cols; j++)
	    {
		FWRITE<<(int)Image[i][j]<<" ";
	    }
	    FWRITE<<endl;
	}
	FWRITE.close();
	return 0;
    }

    int load(string File_Name)
    {
	int PValue,rows,cols;
	string MagicNo;
	string OptComment;
	if (File_Name.substr(File_Name.size()-4,4)!=".pgm")
	{
	    return 1;
	}
	ifstream FREAD(File_Name.c_str());
	if (!FREAD)
	{
	    return 2;
	}
	getline(FREAD,MagicNo);
	getline(FREAD,OptComment);
	FREAD>>cols>>rows>>Maximum;
	setRows(rows);
	setCols(cols);
	for (int i=0; i<Rows; i++)
	{
	    for (int j=0; j<Cols; j++)
	    {
		FREAD>>PValue;
		Image[i][j]=PValue;
	    }
	}
	FREAD.close();
	Loaded=true;
	return 0;
    }

    void combineSideBySide(grayImage &Two, int fillValue = 0)
    {
    	if (!Loaded || !Two.Loaded)
    	{
	    cout<<"Both images must be loaded in order to combine!"<<endl;
	    return;
    	}
	for (int i=0; i<Two.Rows; i++)
	{
	    for (int j=0; j<Two.Cols; j++)
	    {
		if (Cols+j<MaxCols)
		{
		    Image[i][j+Cols]=Two.Image[i][j];
		}
	    }
	}
	if (Rows>Two.Rows)
	{
	    int LastI=MaxCols;
	    if (Cols+Two.Cols<MaxCols)
		LastI=MaxCols-Cols;
	    Fill(Cols,Two.Rows,LastI,Rows-1,fillValue);
	}
	else
	    Fill(0,Rows,Cols-1,Two.Rows-1,fillValue);
	if (Rows<Two.Rows)
	{
	    Rows=Two.Rows;
	}
	Cols+=Two.Cols;
	if (Cols>MaxCols)
	{
	    Cols=MaxCols;
	}
	if (Two.Maximum>Maximum)
	{
	    Maximum=Two.Maximum;
	}
    }

    void combineTopToBottom(grayImage &Two, int fillValue = 0)
    {
    	if (!Loaded || !Two.Loaded)
    	{
	    cout<<"Both images must be loaded in order to combine!"<<endl;
	    return;
    	}
	for (int i=0; i<Two.Rows; i++)
	{
	    for (int j=0; j<Two.Cols; j++)
	    {
		if (Rows+i<MaxRows)
		{
		    Image[i+Rows][j]=Two.Image[i][j];
		}
	    }
	}
	if (Cols>Two.Cols)
	{
	    int LastI=MaxRows-1;
	    if (Rows+Two.Rows<MaxRows)
		LastI=MaxRows-Rows;
	    Fill(Two.Cols,Rows,Cols-1,LastI,fillValue);
	}
	else
	    Fill(Cols,0,Two.Cols-1,Rows-1,fillValue);
	if (Cols<Two.Cols)
	{
	    Cols=Two.Cols;
	}
	Rows+=Two.Rows;
	if (Rows>MaxRows)
	{
	    Rows=MaxRows;
	}
	if (Two.Maximum>Maximum)
	{
	    Maximum=Two.Maximum;
	}
    }

    void FadeIn(grayImage &Two, int Seconds, string BaseFileName, int Frames = 30)
    {
    	if (!Loaded || !Two.Loaded)
    	{
	    cout<<"Both images must be loaded in order to fade in!"<<endl;
	    return;
    	}
        grayImage A;
        double steps=1.0/(Frames*Seconds);
        int R=Rows;
        if(Rows<Two.Rows)
            R=Two.Rows;
        int c= Cols;
        if(Cols<Two.Cols)
            c=Two.Cols;
        A.Rows= R;
        A.Cols= c;
        A.Maximum= Maximum;
	A.Loaded=true;
        int counter=0;
        for(double Alpha=1; Alpha>=0;Alpha-=steps)
        {
            for(int i =0; i<R; i++)
            {
                for(int j=0; j<c; j++)
                    A.Image[i][j]= Image[i][j]*Alpha+(1-Alpha)*Two.Image[i][j];
            }
            stringstream ss;
            ss << counter;
            string str = ss.str();
	    string FileName=BaseFileName+str+".pgm";
            A.Save("Results\\" + FileName);
            counter++;
        }
	
    }

    void Rotate(double angle, int aboutx, int abouty, grayImage &RotatedImage)
    {
	if (!Loaded)
	{
	    cout<<"Image not loaded!"<<endl;
	    return;
	}
	int iNew, jNew;
	angle = (PI*angle)/180.0;
	for (int i=0; i<Rows; i++)
	{
	    for (int j=0; j<Cols; j++)
	    {
		iNew=((i-aboutx)*cos(angle))-((j-abouty)*sin(angle))+aboutx;
		jNew=((i-aboutx)*sin(angle))+((j-abouty)*cos(angle))+abouty;
		if (iNew>0 && iNew<MaxRows && jNew>0 && jNew<MaxCols)
		    RotatedImage.Image[iNew][jNew]=(Image[i][j]+Image[i+1][j]+Image[i][j+1]+Image[i+1][j+1])/4;
	    }
	}
	RotatedImage.Rows=Rows;
	RotatedImage.Cols=Cols;
	RotatedImage.Maximum=Maximum;
	RotatedImage.Loaded=true;
    }

    void Flip(int HorizontalOrVertical=0)
    {
    	if (!Loaded)
    	{
	    cout<<"The image must be loaded in order to flip!"<<endl;
	    return;
    	}
	if (HorizontalOrVertical==0)
	{
	    FlipHorizontal();
	}
	else if (HorizontalOrVertical==1)
	{
	    FlipVertical();
	}
    }

    void Negative()
    {
	if (!Loaded)
	{
	    cout<<"Image not loaded!"<<endl;
	    return;
	}
	for (int i=0; i<Rows; i++)
    	{
	    for (int j=0; j<Cols; j++)
	    {
	    	Image[i][j]=Maximum-Image[i][j];
	    }
    	}
    }

    void changeBrightness(int amount)
    {
	if (!Loaded)
	{
	    cout<<"Image not loaded!"<<endl;
	    return;
	}
	int NewValue;
	for (int i=0; i<Rows; i++)
    	{
	    for (int j=0; j<Cols; j++)
	    {
	    	NewValue=Image[i][j]+amount;
		if (NewValue>Maximum)
		    NewValue=Maximum;
		if (NewValue<0)
		    NewValue=0;
		Image[i][j]=NewValue;
	    }
    	}
    }

    void Quantize(int Levels)
    {
	if (!Loaded)
	{
	    cout<<"Image not loaded!"<<endl;
	    return;
	}
	int Alpha,Step=Maximum/Levels;
	for (int x=0; x<Maximum-Step; x+=Step)
	{
	    Alpha=x+Step;
	    for(int i =0; i<Rows; i++)
            {
            	for(int j=0; j<Cols; j++)
		{
		    if (Image[i][j]>=x && Image[i][j]<Alpha)
		    {
                    	Image[i][j]= (x+Alpha)/2;
		    }
		}
            }
	}
    }

    void medianFilter(grayImage &Result, int filterSize = 3)
    {
	if (!Loaded)
	{
	    cout<<"Image not loaded!"<<endl;
	    return;
	}
	if (filterSize<3 || filterSize%2==0)
	{
	    cout<<"Invalid filter size!";
	    return;
	}
	int A[MaxRows], count;
	for (int r=filterSize; r<(Rows-filterSize); r++)
	{
	    for (int c=filterSize; c<(Cols-filterSize); c++)
	    {
		count=0;
		for (int i=r-filterSize/2; i<=r+filterSize/2; i++)
		{
		    for (int j=c-filterSize/2; j<=c+filterSize/2; j++)
		    {
			A[count]=Image[i][j];
			count++;
		    }
		}
		Sort_Array(A, count);
		Result.Image[r][c]=A[count/2];
	    }
	}
        Result.Rows=Rows;
    	Result.Cols=Cols;
    	Result.Maximum=Maximum;
    	Result.Loaded=true;
    	return;
    }

    void meanFilter(grayImage &Result, int filterSize = 3)
    {
	if (!Loaded)
	{
	    cout<<"Image not loaded!"<<endl;
	    return;
	}
	if (filterSize<3 || filterSize%2==0)
	{
	    cout<<"Invalid filter size!";
	    return;
	}
	int A[MaxRows], count;
	for (int r=filterSize; r<(Rows-filterSize); r++)
	{
	    for (int c=filterSize; c<(Cols-filterSize); c++)
	    {
		count=0;
		for (int i=r-filterSize/2; i<=r+filterSize/2; i++)
		{
		    for (int j=c-filterSize/2; j<=c+filterSize/2; j++)
		    {
			A[count]=Image[i][j];
			count++;
		    }
		}
		int sum, Mean;
		sum=0;
		for (int i=0; i<count; i++)
		{
		    sum+=A[i];
		}
		Mean=sum/count;
		Result.Image[r][c]=Mean;
	    }
	}
    	Result.Rows=Rows;
    	Result.Cols=Cols;
    	Result.Maximum=Maximum;
    	Result.Loaded=true;
    	return;
    }

    void Resize(grayImage &Result,int NewRows, int NewColumns)
    {
	if (!Loaded)
	{
	    cout<<"Image not loaded!"<<endl;
	    return;
	}
        double x=(double)Rows/NewRows;
        double y=(double)Cols/NewColumns;

        for(double i=0,in=0; i<Rows; i+=x,in++)
        {
            for(double j=0,jn=0; j<Cols;j+=y,jn++){
            int I=i, a=in;
            int J=j, b=jn;
            Result.Image[a][b]=Image[I][J];
        }
	}
        Result.Rows= NewRows;
        Result.Cols=NewColumns;
        Result.Maximum= Maximum;
	Result.Loaded=true;
    }

    void Resize(double Ratio,grayImage &result)
    {
	if (!Loaded)
	{
	    cout<<"Image not loaded!"<<endl;
	    return;
	}
        double NewRows=Rows*Ratio;
        double NewColumns = Cols*Ratio;
	if (NewRows>MaxRows)
	    NewRows=MaxRows;
	if (NewColumns>MaxCols)
	    NewColumns=MaxCols;
	int R=NewRows,C=NewColumns;
	Resize(result,R,C);
    }

    void Transform(grayImage &Result, double Matrix[3][3])
    {
	if (!Loaded)
	{
	    cout<<"Image not loaded!"<<endl;
	    return;
	}
        for (int i=0; i<Rows; i++)
	{
	    for (int j=0; j<Cols; j++)
	    {
		double dI=i,dJ=j;
		double I=(dI*Matrix[0][0]) + (dJ*Matrix[0][1]) + Matrix[0][2];
		double J=(dI*Matrix[1][0]) + (dJ*Matrix[1][1]) + Matrix[1][2];
		double Z=(dI*Matrix[2][0]) + (dJ*Matrix[2][1]) + Matrix[2][2];
		if (Z!=0)
		{
		    int NI=I/Z;
		    int NJ=J/Z;
		    if (NI>=0 && NI<MaxRows && NJ>=0 && NJ<MaxCols)
			Result.Image[NI][NJ]=Image[i][j];
		}
	    }
	}
	Result.Rows=Rows;
	Result.Cols=Cols;
	Result.Maximum=Maximum;
	Result.Loaded=true;
    }

    void Filter(grayImage &Result, double Mask[3][3])
    {
	if (!Loaded)
	{
	    cout<<"Image not loaded!"<<endl;
	    return;
	}
	double Val, sum;
	int c=0;
	for (int x=0; x<Rows, c>2; x++)
	{
	    int y=0;
	    c=1;
	    Result.Image[x][y]=Image[x][y];
	    y=Rows-1;
	    c++;
	}
	c=0;
	for (int y=0; y<Cols, c>2; y++)
	{
	    int x=0;
	    c=1;
	    Result.Image[x][y]=Image[x][y];
	    x=Rows-1;
	    c++;
	}
	for (int r=1; r<(Rows-1); r++)
	{
	    for (int c=1; c<(Cols-1); c++)
	    {
		sum=0;
		for (int i=r-1, x=0; i<=r+1; i++, x++)
		{
		    for (int j=c-1,y=0; j<=c+1; j++, y++)
		    {
			Val=(double)Image[i][j];
			sum=sum+Val*Mask[x][y];
		    }
		}
		if (sum>Maximum)
		    sum=Maximum;
		if (sum<0)
		    sum=0;
		int NewVal=sum;
		Result.Image[r][c]=NewVal;
	    }
	}
	Result.Maximum=Maximum;
	Result.Rows=Rows;
	Result.Cols=Cols;
	Result.Loaded=true;
    }

    void DerivativeImage(grayImage &Result)
    {
    	if (!Loaded)
	{
	    cout<<"Image not loaded!"<<endl;
	    return;
	}
	int V1[3][3], V2[3][3];
	int Val1,Val2=-1;
	int Val,sum1,sum2,FinalVal;
	for (int r=0; r<3; r++)
	{
	    Val1=-1;
	    for (int c=0; c<3; c++)
	    {
		V2[r][c]=Val2;
		V1[r][c]=Val1;
		Val1++;
	    }
	    Val2++;
	}
	int c=0;
	for (int x=0; x<Rows, c>2; x++)
	{
	    int y=0;
	    c=1;
	    Result.Image[x][y]=Image[x][y];
	    y=Rows-1;
	    c++;
	}
	c=0;
	for (int y=0; y<Cols, c>2; y++)
	{
	    int x=0;
	    c=1;
	    Result.Image[x][y]=Image[x][y];
	    x=Rows-1;
	    c++;
	}
	for (int r=1; r<(Rows-1); r++)
	{
	    for (int c=1; c<(Cols-1); c++)
	    {
		sum1=0;
		sum2=0;
		for (int i=r-1, x=0; i<=r+1; i++, x++)
		{
		    for (int j=c-1,y=0; j<=c+1; j++, y++)
		    {
			Val=Image[i][j];
			sum1=sum1+Val*V1[x][y];
			sum2=sum2+Val*V2[x][y];
		    }
		}
		FinalVal=sqrt((sum1*sum1)+(sum2*sum2));
		if (FinalVal<0)
		    FinalVal=0;
		if (FinalVal>255)
		    FinalVal=255;
		Result.Image[r][c]=FinalVal;
	    }
	}
	Result.Rows=Rows;
	Result.Cols=Cols;
	Result.Maximum=255;
	Result.Loaded=true;
    }

    void Crop(grayImage &result,int L, int T, int R, int B, int RSF=0)
    {
	if (!Loaded)
	{
	    cout<<"Image not loaded!"<<endl;
	    return;
	}
	if (T<0 || L<0 || B>Rows || R>Cols)
	{
	    cout<<"Invalid coordinates entered!";
	    return;
	}
        int i,j,x,y;
	for( i=0,x=T; x<B;i++,x++)
        {
            for(j=0,y=L; y<R;j++,y++)
            	result.Image[i][j]= Image[x][y];
        }
        result.Rows=(B-T);
        result.Cols=(R-L);
        result.Maximum= Maximum;
        result.Loaded=true;
        if(RSF==1)
	{
	    grayImage Two;
            result.Resize(Two,Rows,Cols);
	    for (int i=0; i<Rows; i++)
		for (int j=0; j<Cols; j++)
		    result.Image[i][j]=Two.Image[i][j];
	    result.Rows=Rows;
	    result.Cols=Cols;
	}
    }

private:

    void FlipHorizontal(){

	int temp;
	for (int i=0; i<Rows/2; i++)
    	{
	    for (int j=0; j<Cols; j++)
	    {
	    	temp=Image[i][j];
		Image[i][j]=Image[Rows-1-i][j];
		Image[Rows-1-i][j]=temp;
	    }
    	}
    }
    void FlipVertical(){

	int temp;
	for (int j=0; j<Cols/2; j++)
    	{
	    for (int i=0; i<Rows; i++)
	    {
	    	temp=Image[i][j];
		Image[i][j]=Image[i][Cols-1-j];
		Image[i][Cols-1-j]=temp;
	    }
    	}
    }

    void Sort_Array(int Arr[], int n)
    {
    	int Previous,temp;
    	for (int j=0; j<n-1; j++)
    	{
            Previous=j;
            for (int i=j+1; i<n; i++)
            {
        	if (Arr[i]<Arr[Previous])
        	{
               	    Previous=i;
        	}
                temp=Arr[j];
        	Arr[j]=Arr[Previous];
        	Arr[Previous]=temp;
    	    }
    	}
        return;
    }

    void Fill(int L, int T, int R, int B, int FillValue)
    {
	if (FillValue>255)
	    FillValue=255;
	if (FillValue<0)
	    FillValue=0;
        for(int i = T; i<= B; i++)
            for(int j = L; j <= R; j++)
                Image[i][j] = FillValue;
    }

    unsigned short Image[MaxRows][MaxCols];
    int Rows, Cols, Maximum;
    bool Loaded;
};

void inputMatrix(double Matrix[3][3])
{
    for (int i=0; i<3; i++)
    {
	for (int j=0; j<3; j++)
	{
	    cin>>Matrix[i][j];
	}
    }
}

int main()
{
    grayImage GM;
    double Arr[3][3];
    string FileName;
    int choice;
    do
    {
	grayImage GM2;
	cout<<"1. Load an image."<<endl;
	cout<<"2. Save the new image."<<endl;
	cout<<"3. Combine two images."<<endl;
	cout<<"4. Fade one image to form another image."<<endl;
	cout<<"5. Rotate an image."<<endl;
	cout<<"6. Flip an image."<<endl;
	cout<<"7. Take negative of image."<<endl;
	cout<<"8. Change the brightness of image."<<endl;
	cout<<"9. Quantize image."<<endl;
	cout<<"10. Apply median filter to image."<<endl;
	cout<<"11. Apply mean filter to image."<<endl;
	cout<<"12. Resize an image."<<endl;
	cout<<"13. Crop an image."<<endl;
	cout<<"14. Apply a filter using a mask on image."<<endl;
	cout<<"15. Take Derivative of image."<<endl;
	cout<<"16. Transform an image using transformation matrix."<<endl;
	cout<<"Enter -999 to terminate."<<endl<<endl;;
	cout<<"Enter choice[1..16]: ";
	cin>>choice;
	if (choice==1)
	{
	    cout<<"Enter name of file to load: ";
	    cin>>FileName;
	    int flag=GM.load(FileName);
	    if (flag==1)
		cout<<"Image extension must be .pgm!"<<endl;
	    else if (flag==2)
		cout<<"File could not be opened!"<<endl;
	    else
		cout<<"Image loaded successfully!"<<endl;
	}
	else if (choice==2)
	{
	    cout<<"Enter name of new file to save as: ";
	    cin>>FileName;
	    int flag=GM.Save(FileName);
	    if (flag==1)
		cout<<"Image extension must be .pgm!"<<endl;
	    else if (flag==2)
		cout<<"File could not be opened!"<<endl;
	    else if (flag==-1)
	        cout<<"Image needs to be loaded first!"<<endl;
	    else
		cout<<"New image saved successfully!"<<endl;
	}
        else if (choice==3)
	{
	    int choice2,FVal;
	    cout<<"1. Combine side by side."<<endl;
	    cout<<"2. Combine top to bottom."<<endl; 
	    cout<<"Enter choice[1-2]: ";
	    cin>>choice;
	    if (choice<1 || choice>2)
	    {
		cout<<"Invalid choice entered!"<<endl;
	    }
	    else
	    {
		cout<<"Enter name of second image to load: ";
		cin>>FileName;
		int flag1=GM2.load(FileName);
		if (flag1==1)
		    cout<<"Image extension must be .pgm!"<<endl;
	    	else if (flag1==2)
		    cout<<"File could not be opened!"<<endl;
		else
		{
	    	    cout<<"Enter fill value: ";
	    	    cin>>FVal;
		    cout<<"Enter name of combined image to save: ";
		    cin>>FileName;
	    	    if (choice==1)
		    	GM.combineSideBySide(GM2,FVal);
	    	    else if (choice==2)
		    	GM.combineTopToBottom(GM2,FVal);
		    int flag2=GM.Save(FileName);
	    	    if (flag2==1)
	  	    	cout<<"Image extension must be .pgm!"<<endl;
	    	    else if (flag2==2)
		    	cout<<"File could not be opened!"<<endl;
	    	    else if (flag2==0)
		    	cout<<"Image saved successfully!"<<endl;
		}
	    }
	}
        else if (choice==4)
	{
	    int Seconds, Frames;
	    cout<<"Enter name of second image to load: ";
	    cin>>FileName;
	    int flag1=GM2.load(FileName);
	    if (flag1==1)
		cout<<"Image extension must be .pgm!"<<endl;
	    else if (flag1==2)
		cout<<"File could not be opened!"<<endl;
	    else
	    {
	    	cout<<"Enter seconds: ";
	    	cin>>Seconds;
	    	cout<<"Enter number of frames: ";
	    	cin>>Frames;
	    	cout<<"Enter base file name to save as (without extension): ";
	    	cin>>FileName;
	    	GM.FadeIn(GM2, Seconds, FileName, Frames);
		cout<<"Images created successfully!"<<endl;
	    }
	}
	else if (choice==5)
	{
	    int DegAngle,AboutX,AboutY;
	    cout<<"Enter angle in degrees to rotate: ";
	    cin>>DegAngle;
	    cout<<"Enter coordinates to rotate about: "<<endl;
	    cout<<"Enter x-coordinate: ";
	    cin>>AboutX;
	    cout<<"Enter y-coordinate: ";
	    cin>>AboutY;
	    GM.Rotate(DegAngle, AboutX, AboutY, GM2);
	    cout<<"Enter name of file to save rotated image as: ";
	    cin>>FileName;
	    int flag=GM2.Save(FileName);
	    if (flag==1)
		cout<<"Image extension must be .pgm!"<<endl;
	    else if (flag==2)
		cout<<"File could not be opened!"<<endl;
	    else if (flag==0)
		cout<<"Rotated image saved successfully!"<<endl;
	}
	else if (choice==6)
	{
	    int choice2;
	    cout<<"0. Flip Horizontally."<<endl;
	    cout<<"1. Flip Vertically."<<endl;
	    cout<<"Enter choice[0-1]: ";
	    cin>>choice2;
	    if (choice2==0 || choice2==1)
	    {
		cout<<"Enter name of file to save flipped image: ";
		cin>>FileName;
	    	GM.Flip(choice2);
		int flag2=GM.Save(FileName);
	    	if (flag2==1)
	  	    cout<<"Image extension must be .pgm!"<<endl;
	    	else if (flag2==2)
		    cout<<"File could not be opened!"<<endl;
	    	else if (flag2==0)
		    cout<<"Image saved successfully!"<<endl;
	    }
	    else
		cout<<"Invalid choice entered!"<<endl;
	}
	else if (choice==7)
	{
	    cout<<"Enter name of file to save negative image: ";
	    cin>>FileName;
	    GM.Negative();
	    int flag2=GM.Save(FileName);
	    if (flag2==1)
	  	cout<<"Image extension must be .pgm!"<<endl;
	    else if (flag2==2)
		cout<<"File could not be opened!"<<endl;
	    else if (flag2==0)
		cout<<"Image saved successfully!"<<endl;
	}
	else if (choice==8)
	{
	    int amt;
	    cout<<"Enter amount to change brightness: ";
	    cin>>amt;
	    cout<<"Enter name of file to save brightened/darkened image: ";
	    cin>>FileName;
	    GM.changeBrightness(amt);
	    int flag2=GM.Save(FileName);
	    if (flag2==1)
	  	cout<<"Image extension must be .pgm!"<<endl;
	    else if (flag2==2)
		cout<<"File could not be opened!"<<endl;
	    else if (flag2==0)
		cout<<"Image saved successfully!"<<endl;
	}
	else if (choice==9)
	{
	    int lvl;
	    cout<<"Enter number of levels to quantize: ";
	    cin>>lvl;
	    cout<<"Enter name of file to save quantized image: ";
	    cin>>FileName;
	    GM.Quantize(lvl);
	    int flag2=GM.Save(FileName);
	    if (flag2==1)
	  	cout<<"Image extension must be .pgm!"<<endl;
	    else if (flag2==2)
		cout<<"File could not be opened!"<<endl;
	    else if (flag2==0)
		cout<<"Image saved successfully!"<<endl;
	}
	else if (choice==10)
	{
	    int fsize;
	    cout<<"Enter filter size: ";
	    cin>>fsize;
	    cout<<"Enter name of file to save modified image: ";
	    cin>>FileName;
	    GM.medianFilter(GM2,fsize);
	    int flag2=GM2.Save(FileName);
	    if (flag2==1)
	  	cout<<"Image extension must be .pgm!"<<endl;
	    else if (flag2==2)
		cout<<"File could not be opened!"<<endl;
	    else if (flag2==0)
		cout<<"Image saved successfully!"<<endl;
	}
	else if (choice==11)
	{
	    int fsize;
	    cout<<"Enter filter size: ";
	    cin>>fsize;
	    cout<<"Enter name of file to save modified image: ";
	    cin>>FileName;
	    GM.meanFilter(GM2,fsize);
	    int flag2=GM2.Save(FileName);
	    if (flag2==1)
	  	cout<<"Image extension must be .pgm!"<<endl;
	    else if (flag2==2)
		cout<<"File could not be opened!"<<endl;
	    else if (flag2==0)
		cout<<"Image saved successfully!"<<endl;
	}
	else if (choice==12)
	{
	    int choice2, NRows,NCols;
	    double Ratio;
	    cout<<"1. Resize using ratio."<<endl;
	    cout<<"2. Resize by entering size of new rows and new columns."<<endl;
	    cout<<"Enter choice[1-2]: ";
	    cin>>choice2;
	    if (choice2==1 || choice2==2)
	    {
		if (choice2==1)
		{
		    cout<<"Enter ratio: ";
		    cin>>Ratio;
		    cout<<"Enter name of resized image to save as: ";
		    cin>>FileName;
		    GM.Resize(Ratio, GM2);
		}
		else
		{
		    cout<<"Enter new rows: ";
		    cin>>NRows;
		    cout<<"Enter new columns: ";
		    cin>>NCols;
		    cout<<"Enter name of resized image to save as: ";
		    cin>>FileName;
		    GM.Resize(GM2, NRows,NCols);
		}
	    	int flag2=GM2.Save(FileName);
	    	if (flag2==1)
	  	    cout<<"Image extension must be .pgm!"<<endl;
	    	else if (flag2==2)
		    cout<<"File could not be opened!"<<endl;
	    	else if (flag2==0)
		    cout<<"Image saved successfully!"<<endl;
	    }
	    else
	    {
		cout<<"Invalid choice entered!"<<endl;
	    }
	}
	else if (choice==13)
	{
	    int Top,Bottom,Right,Left,RSF;
	    cout<<"Enter coordinates of image from where it is to be cropped: "<<endl;
	    cout<<"Enter top coordinate: ";
	    cin>>Top;
	    cout<<"Enter bottom coordinate: ";
	    cin>>Bottom;
	    cout<<"Enter left coordinate: ";
	    cin>>Left;
	    cout<<"Enter top coordinate: ";
	    cin>>Right;
	    cout<<"Enter 1 if the image should be resized to its original size and 0 if not: ";
	    cin>>RSF;
	    cout<<"Enter name of file to save cropped image as: ";
	    cin>>FileName;
	    GM.Crop(GM2,Left,Top,Right,Bottom,RSF);
	    int flag2=GM2.Save(FileName);
	    if (flag2==1)
	  	cout<<"Image extension must be .pgm!"<<endl;
	    else if (flag2==2)
		cout<<"File could not be opened!"<<endl;
	    else if (flag2==0)
		cout<<"Image saved successfully!"<<endl;
	}
	else if (choice==14)
	{
	    cout<<"Input values of mask: ";
	    inputMatrix(Arr);
	    cout<<"Enter name of file to save modified image as: ";
	    cin>>FileName;
	    GM.Filter(GM2, Arr);
	    int flag2=GM2.Save(FileName);
	    if (flag2==1)
	  	cout<<"Image extension must be .pgm!"<<endl;
	    else if (flag2==2)
		cout<<"File could not be opened!"<<endl;
	    else if (flag2==0)
		cout<<"Image saved successfully!"<<endl;
	}
	else if (choice==15)
	{
	    cout<<"Enter name of file to save derivative image as: ";
	    cin>>FileName;
	    GM.DerivativeImage(GM2);
	    int flag2=GM2.Save(FileName);
	    if (flag2==1)
	  	cout<<"Image extension must be .pgm!"<<endl;
	    else if (flag2==2)
		cout<<"File could not be opened!"<<endl;
	    else if (flag2==0)
		cout<<"Image saved successfully!"<<endl;
	}
	else if (choice==16)
	{
	    cout<<"Input values of transformation matrix: ";
	    inputMatrix(Arr);
	    cout<<"Enter name of file to save transformed image as: ";
	    cin>>FileName;
	    GM.Transform(GM2, Arr);
	    int flag2=GM2.Save(FileName);
	    if (flag2==1)
	  	cout<<"Image extension must be .pgm!"<<endl;
	    else if (flag2==2)
		cout<<"File could not be opened!"<<endl;
	    else if (flag2==0)
		cout<<"Image saved successfully!"<<endl;
	}
	else if (choice!=-999)
	{
	    cout<<"Invalid choice entered!"<<endl;
	}
	cout<<endl<<endl;
    }
    while (choice!=-999);

    return 0;
}