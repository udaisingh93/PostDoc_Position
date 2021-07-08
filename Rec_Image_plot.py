##################################################
# Author: {Udai Singh}
## Email: {udai.singh93@gmail.com}
## Status: {development}
##################################################
import os
##################################################
#small function to install  and check libraries
##################################################
def install_or_import_package(package):
    try:
        abv= __import__(package)
    except:
        print(f'{package} not found on system trying to install it')
        os.system(f'pip3 install {package}')
        abv= __import__(package)
    return abv
##################################################
#importing required libraries
##################################################
pd=install_or_import_package('pandas')
np=install_or_import_package('numpy')
matplotlib=install_or_import_package('matplotlib')
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
argparse=install_or_import_package('argparse')
##################################################
#defining variables
##################################################
inputfile=""
outputfile="RecImageCor.dat"
sepr="\s+"
skiprow=0
##################################################
#parsing variables
##################################################
parser=argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str,help="Enter input file name")
parser.add_argument("-o", "--output", type=str,help="Enter output file name")
parser.add_argument("-s", "--sep", type=str,help="Enter column seprater as str")
parser.add_argument("-sl", "--skiplines", type=int,help="Enter no. of lines to be skiped")
args = parser.parse_args()
if(args.input==None):
    print("Please enter input flag and file name")
    print("Please enter flag -h for help")
    exit()
else:
    inputfile=args.input
    
if(args.output!=None):
    outputfile=args.output
if(args.sep!=None):
    sepr=args.sep
    print(f'Seperator changed to {sepr}')
if(args.skiplines!=None):
    skiprow=args.skiplines
    print(f'Skiping {skiprow} lines from input file')
print(f'Input file is {inputfile}')
print(f'Output file is {outputfile}')
####################################################################################################
#Reading input file and  coverting it data frame and numpy arrays
####################################################################################################
df=pd.read_csv(inputfile,header=None,sep=sepr,skiprows=skiprow) #read ASCII file convert to pandas data frame
#list of parameter read from ascii file will be assigned to dataframe
hedList=['globalPosX1','globalPosY1', 'globalPosZ1','time1','globalPosX2','globalPosY2', 'globalPosZ2','time2']
ListAnnihilaton=['X','Y','Z'] #List of parameter for output
pointA=np.array(df.loc[:,0:2]) #point1 one numpy array
pointB=np.array(df.loc[:,4:6]) #point2 two nu
####################################################################################################
#Calculating annihilation point
####################################################################################################
df.columns=hedList #assign columns name to data frame
df['shift']=(df['time1']-df['time2'])*299.792*0.5 #calculate shift 
midleOfLor=(pointA+pointB)/2.0
versorOnLor=(pointA-pointB)
annihilationPoint=midleOfLor+np.array(df['shift']).reshape(-1,1)*versorOnLor
df2=pd.DataFrame(data=annihilationPoint,columns=ListAnnihilaton)
df.to_csv(outputfile,index=None,sep='\t') #saving data to output file
####################################################################################################
#ploting and saving figure 
####################################################################################################
fontsize=20
plt.figure(figsize=(10,10))
plt.hist2d(df2['X'],df2['Y'],range=[[-15,15],[-15,15]],bins=50,cmap='plasma',norm=LogNorm())
plt.title("Reconstructed source image",fontsize=fontsize)
plt.xlabel("X [mm]",fontsize=fontsize)
plt.ylabel("Y [mm]",fontsize=fontsize)
plt.savefig("sourceimage.png")