#/bin/sh

### off-centering model: [gauss,nfw]
offmodel=gauss

### unit of spatial coordinates [comoving, phys]
unit=phys
#unit=comoving

### directory of measurement data
datadir=data/

### number of MCMC samples
nsample=4000000

### add the prior of minimum qcen value for high-pcen sample? 
qcen_min_highpcen=0.95
#qcen_min_highpcen=0

### sample [boss,dr7,new]
sample=new

### measurement [angcor,dsigma]
measure=angcor
#measure=dsigma

### covariance is included ? [y,n]
ccov=y

#!!!
if [ $measure = 'dsigma' ]; then
ccov=n
fi


###ndata='single double triple'
ndata=triple

### comoving scale range for fitting [unit: Mpc/h]
rmin=0.1
rmax=2

### use generalized NFW profile ? [y,n]
igenNFW=n

### 2halo term is added? [y,n]
i2halo=n

### output filename
outf=chains/${measure}_${sample}_rmax${rmax}_rmin${rmin}_${ndata}_${unit}

if [ ! -d chains/ ]; then
mkdir chains
fi

if [ ! -d mpiout/ ]; then
mkdir mpiout
fi

if [ $igenNFW = 'y' ]; then
outf=${outf}_genNFW
fi

if [ $i2halo = 'y' ]; then
outf=$outf'_w2h'
fi

if [ $ccov = 'n' ]; then
outf=$outf'_nocov'
fi

outf=${outf}_${offmodel}

echo 'output filename: '$outf

#######################################################################################

file='modelfit_in.f90'
filesub='modelfit.f90'

echo 'MODULE Readmodel' > $file
echo '  implicit none' >> $file
echo "  character(len=6), parameter :: measurement='"${measure}"'" >> $file
echo '  real :: rmin='$rmin', rmax='$rmax >> $file

if [ $ndata = 'single' ]; then
echo "  integer, parameter :: nlrg=1" >> $file
echo "  character(len=10), parameter :: clrg(1)=(/'highpcen'/)" >> $file
echo "  character(len=100), parameter :: infname='"${datadir}"redmapper."${measure}"."${unit}"."${sample}".'" >> $file
elif [ $ndata = 'double' ]; then
echo "  integer, parameter :: nlrg=2" >> $file
echo "  character(len=10), parameter :: clrg(2)=(/'rmcg','blrg'/)" >> $file
echo "  character(len=100), parameter :: infname='"${datadir}"redmapper."${measure}"."${unit}"."${sample}".'" >> $file
elif [ $ndata = 'triple' ]; then
echo "  integer, parameter :: nlrg=3" >> $file
echo "  character(len=10), parameter :: clrg(3)=(/'highpcen','rmcg','blrg'/)" >> $file
echo "  character(len=100), parameter :: infname='"${datadir}"redmapper."${measure}"."${unit}"."${sample}".'" >> $file
fi

echo "  character(len=100), parameter :: inf2h='"${datadir}"2halo_"${measure}".dat'" >> $file
echo "  character(len=1), parameter :: ccov='"${ccov}"',genNFW='"${igenNFW}"'" >> $file
echo "  character(len=5), parameter :: offmodel='"${offmodel}"'" >> $file
echo "  double precision, parameter :: pi=3.1415926535897932384" >> $file
if [ $unit = 'phys' ]; then
echo "  double precision, parameter :: z=0.25,omegam=0.3,gconst=4.3e-9,rhoc=3e4/(8*pi*gconst),rhom=rhoc*omegam*(1+z)**3">> $file
else
echo "  double precision, parameter :: z=0.25,omegam=0.3,gconst=4.3e-9,rhoc=3e4/(8*pi*gconst),rhom=rhoc*omegam" >> $file
fi

cat $filesub >> $file

make clean
make cosmomc

#######################################################################################

paramfile='params.ini'
paramfile_sub='params_sub.ini'

echo 'file_root = '$outf > $paramfile
echo 'samples = '$nsample >> $paramfile
cat $paramfile_sub >> $paramfile

################################ parameter range #####################################

if [ $measure = 'dsigma' ]; then

#param1: Average host halo mass M180b [unit:10^14Msun/h]
echo 'param1 = 2. 0.1 10. 0.01 0.01' >> $paramfile

#param2: concentration parameter
echo 'param2 = 5. 0.01 20. 0.01 0.01' >> $paramfile

else

#param1: Amplitude of 1-halo term
echo 'param1 = 7 0 50 0.01 0.01' >> $paramfile
#echo 'param1 = 1 0 10 0.01 0.01' >> $paramfile

#param2: r_s: transition scale of NFW profile
#echo 'param2 = 0.2 0 5 0.001 0.001' >> $paramfile

#param2: c_amp: transition scale of NFW profile
echo 'param2 = 1 0 5 0.01 0.01' >> $paramfile

fi

#param3-5: offcentering scale for 3 samples [unit:Mpc/h] 
if [ $ndata = 'single' ]; then
echo 'param3 = 0.4 0. 1. 0.002 0.002' >> $paramfile
echo 'param4 = 0 0 0 0 0' >> $paramfile
echo 'param5 = 0 0 0 0 0' >> $paramfile

elif [ $ndata = 'double' ]; then
echo 'param3 = 0.4 0. 1. 0.002 0.002' >> $paramfile
echo 'param4 = 0.4 0. 1. 0.002 0.002' >> $paramfile
echo 'param5 = 0 0 0 0 0' >> $paramfile

elif [ $ndata = 'triple' ]; then
if [ $offmodel = 'nfw' ]; then
echo 'param3 = 1 0 10 0.3 0.3' >> $paramfile
echo 'param4 = 1 0 10 0.3 0.3' >> $paramfile
echo 'param5 = 1 0 10 0.3 0.3' >> $paramfile
else
echo 'param3 = 0.4 0. 1. 0.002 0.002' >> $paramfile
echo 'param4 = 0.4 0. 1. 0.002 0.002' >> $paramfile
echo 'param5 = 0.4 0. 1. 0.002 0.002' >> $paramfile
fi
fi

#param6-8: subhalo mass for 3 samples [unit:10^14Msun/h]
if [ $measure = 'dsigma' ]; then
if [ $ndata = 'single' ]; then
echo 'param6 = 0.02 0 0.1 0.0002 0.0002' >> $paramfile
echo 'param7 = 0 0 0 0 0' >> $paramfile
echo 'param8 = 0 0 0 0 0' >> $paramfile
elif [ $ndata = 'double' ]; then
echo 'param6 = 0.02 0 0.1 0.0002 0.0002' >> $paramfile
echo 'param7 = 0.02 0 0.1 0.0002 0.0002' >> $paramfile
echo 'param8 = 0 0 0 0 0' >> $paramfile
elif [ $ndata = 'triple' ]; then
echo 'param6 = 0.02 0 0.1 0.0002 0.0002' >> $paramfile
echo 'param7 = 0.02 0 0.1 0.0002 0.0002' >> $paramfile
echo 'param8 = 0.02 0 0.1 0.0002 0.0002' >> $paramfile
fi
else
echo 'param6 = 0 0 0 0 0' >> $paramfile
echo 'param7 = 0 0 0 0 0' >> $paramfile
echo 'param8 = 0 0 0 0 0' >> $paramfile
fi

#param9-11: central fraction for 3 samples
echo 'param9 = 0.97 '$qcen_min_highpcen' 1 0.002 0.002' >> $paramfile

if [ $ndata = 'single' ]; then
echo 'param10 = 0 0 0 0 0' >> $paramfile
echo 'param11 = 0 0 0 0 0' >> $paramfile
elif [ $ndata = 'double' ]; then
echo 'param10 = 0.3 0 1 0.002 0.002' >> $paramfile
echo 'param11 = 0 0 0 0 0' >> $paramfile
elif [ $ndata = 'triple' ]; then
echo 'param10 = 0.75 0 1 0.002 0.002' >> $paramfile
echo 'param11 = 0.16 0 1 0.002 0.002' >> $paramfile
fi

#param12: amplitude of 2-halo term
if [ $i2halo = 'y' ]; then
echo 'param12 = 0.2 0.01 2 0.001 0.001' >> $paramfile
else
echo 'param12 = 0 0 0 0 0' >> $paramfile
fi

#param13-14: outer/inner power-law index of generalized NFW profile
if [ $igenNFW = 'y' ]; then
echo 'param13 = 3 1 5 0.01 0.01' >> $paramfile
echo 'param14 = 1 0 3 0.01 0.01' >> $paramfile
else
echo 'param13 = 3 3 3 0 0' >> $paramfile
echo 'param14 = 1 1 1 0 0' >> $paramfile
fi

#param15
echo 'param15 = 0 0 0 0 0' >> $paramfile


