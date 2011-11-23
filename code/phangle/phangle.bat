mkdir phangle$1_$2
cd phangle$1_$2
cp ../rsave.na8-init rsave.na8-phangle

echo "&GLOBAL" > for005.na8-phangle
echo "nclust=8,nion=8," >> for005.na8-phangle
echo "nspdw=4," >> for005.na8-phangle
echo "nion2=0," >> for005.na8-phangle
echo "temp=0.00," >> for005.na8-phangle
echo "dx=0.8,dy=0.8,dz=0.8," >> for005.na8-phangle
echo "kxbox=48,kybox=48,kzbox=48,kstate=12," >> for005.na8-phangle
echo "b2occ=0.2,deocc=0.4," >> for005.na8-phangle
echo "init_lcao=0," >> for005.na8-phangle
echo "epsoro=1.0e-6," >> for005.na8-phangle
echo "epswf=0.1,e0dmp=1.2," >> for005.na8-phangle
echo "&END" >> for005.na8-phangle
echo "&DYNAMIC" >> for005.na8-phangle
echo "iflocaliz=0,jelf=0," >> for005.na8-phangle
echo "ifhamdiag=0," >> for005.na8-phangle
echo "dt1=0.1," >> for005.na8-phangle
echo "istinf=10," >> for005.na8-phangle
echo "jstinf=10,jdip=10," >> for005.na8-phangle
echo "jenergy=10," >> for005.na8-phangle
echo "isave=00,istat=1,irest=0," >> for005.na8-phangle
echo "ismax=00," >> for005.na8-phangle
echo "itmax=10," >> for005.na8-phangle
echo "ifsicp=2," >> for005.na8-phangle
echo "ipsptyp=0," >> for005.na8-phangle
echo "jinfo=10," >> for005.na8-phangle
echo "jquad=10," >> for005.na8-phangle
echo "ifredmas=1," >> for005.na8-phangle
echo "ionmdtyp=0," >> for005.na8-phangle
echo "phangle=$1,phphase=$2,npstate=9,nhstate=4," >> for005.na8-phangle
echo "&END" >> for005.na8-phangle
echo "&SURFACE" >> for005.na8-phangle
echo "iPotStatic=0," >> for005.na8-phangle
echo "&END" >> for005.na8-phangle

echo "na8-phangle" > for005

~/triax/master/essai.seq > for007

echo "$1 $2" >> ../collect
tail -n 1 penergies.na8-phangle >> ../collect

cd ..
