program Hess2freq
implicit real*8 (a-h,o-z)
real*8,parameter :: pi=3.141592653589793D0
real*8,allocatable :: coord(:),atomid(:)
real*8,allocatable :: Hess(:,:),eigvecmat(:,:),eigvalarr(:),kmat(:,:),massmat(:,:),atmmass(:),freq(:),wavenum(:),normvec(:,:),tmpvec(:)
character*200 filename
real*8,allocatable :: Rsamp(:),Psamp(:)
real*8,allocatable :: cov_R(:,:),cov_P(:,:),mean_R(:),mean_P(:),vec(:),mat(:,:)
real*8 :: beta,T
integer :: i_typesamp, i_typeout
logical alive
character*10 c10
character*300 c300
character(len=2), parameter :: atom_names(0:109) = [ &
    "Bq","H ","He",  & ! 0–2
    "Li","Be","B ","C ","N ","O ","F ","Ne",  & ! 3–10
    "Na","Mg","Al","Si","P ","S ","Cl","Ar",  & ! 11–18
    "K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn", &
    "Ga","Ge","As","Se","Br","Kr",  & ! 19–36
    "Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd", &
    "In","Sn","Sb","Te","I ","Xe",  & ! 37–54
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy", &
    "Ho","Er","Tm","Yb","Lu",  & ! 55–71
    "Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi", &
    "Po","At","Rn",  & ! 72–86
    "Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm","Bk","Cf", &
    "Es","Fm","Md","No","Lr",  & ! 87–103
    "Rf","Db","Sg","Bh","Hs","Mt" ]  ! 104–109


!Don't support linear molecules
write(*,*) "Hess2freq: Load Hessian from .fch file and then compute harmonic frequencies"
write(*,*) "Programmed by Sobereva (Sobereva@sina.com)"
write(*,*) "Release date: 2016-May-27"
write(*,*)

isilent=0
call getarg(1,filename)
if (trim(filename)=="") then
	write(*,*) "Input path of a .fch file containing Hessian matrix (e.g. freq task)"
	do while(.true.)
		read(*,"(a)") filename
		inquire(file=filename,exist=alive)
		if (alive) exit
		write(*,*) "Cannot find the file, input again"
	end do
else
	isilent=1
end if

open(10,file=filename,status="old")
call loclabel(10,"Number of atoms")
read(10,"(49x,i12)") natm
nmode=3*natm-6 !The number of vibrational modes
nmodeall=3*natm
write(*,"(' The number of atoms:',i6)") natm
write(*,"(' The number of vibrational modes:',i6)") nmode

!Load atomic ID
allocate(atomid(natm))
call loclabel(10,"Atomic numbers")
read(10,"(a)") c300
! print *,c300
if(natm < 6)then
	read(10,*) atomid
else
	! read(10,"(6(1PI))") atomid
	nline = int(natm / 6) + 1
	nres = mod(natm, 6)
	do i=1,nline
		read(10,*) atomid((i-1)*6+1:min(i*6,natm))
	end do
	if (nres > 0) then
		read(10,*) atomid(nline*6+1:natm)
	end if
end if
write(*,*) "atomic id:"
write(*,"(6(I4))") int(atomid) !Convert to integer
! stop
write(*,*) "atomic list:"
write(*,"(10(a))") atom_names(int(atomid))

! stop


!Load atomic coordinates
allocate(coord(natm*3))
call loclabel(10,"Current cartesian coordinates")
read(10,*)
read(10,"(5(1PE16.8))") coord
write(*,*) "atomic coordinate:"
write(*,"(6(E18.8))") coord * 0.529177249 !Convert to Angstrom
! stop

!Load atomic masses
allocate(atmmass(natm))
call loclabel(10,"Real atomic weights")
read(10,*)
read(10,"(5(1PE16.8))") atmmass
write(*,*) "Atomic masses:"
write(*,"(6(f12.6))") atmmass

atmmass=atmmass*1.82289E3

!Load Hessian matrix
allocate(Hess(nmodeall,nmodeall))
Hess=0
call loclabel(10,"Cartesian Force Constants")
read(10,*)
read(10,"(5(1PE16.8))") ((Hess(i,j),j=1,i),i=1,nmodeall)
Hess=Hess+transpose(Hess)
do i=1,nmodeall
	Hess(i,i)=Hess(i,i)/2D0
end do
close(10)
! call showmatgau(Hess,"Hessian matrix",1)

!Construct mass matrix
allocate(massmat(nmodeall,nmodeall))
massmat=0
do i=1,natm
	do j=(i-1)*3+1,i*3
		massmat(j,j)=atmmass(i)
	end do
end do



!Construct force constant matrix (mass-weighted Hessian matrix)
allocate(kmat(nmodeall,nmodeall))
do i=1,nmodeall
	do j=1,nmodeall
		kmat(i,j)=hess(i,j)/dsqrt(massmat(i,i)*massmat(j,j))
	end do
end do
! write(*,*)
call showmatgau(kmat,"Force constant matrix (i.e mass-weighted Hessian)")

!Diagonalization of force constant matrix
allocate(eigvecmat(nmodeall,nmodeall),eigvalarr(nmodeall))
call diagsymat(kmat,eigvecmat,eigvalarr,istat) !Note the column of eigvecmat has already been automatically normalized
if (istat==0) write(*,*) "Diagonalization passed"

write(*,*) "frequency (cm-1) of each mode:"
write(*,*) sqrt(abs(eigvalarr))*sign(1.0d0,eigvalarr) * 219474.6313702



write(*,"(a)") "Input T (K), Nsamp, type of sampling (0 for classical, 1 for quantum), and output form (0 for 1 xyz file, 1 for multiple xyz files):"
read(*,*) T,nssmp,i_typesamp,i_typeout



beta=1.0D0/(T*1.3806E-5/4.35974) !Beta in a.u.

write(*,*) "T=",T, "K, beta=", beta, ", Nsamp=",nssmp," i_typesamp=",i_typesamp," i_typeout=",i_typeout

if(i_typeout == 0)then 
	open(11,file="Rsample.xyz",status="replace")
	open(12,file="Psample.xyz",status="replace")
end if
allocate(Rsamp(3*natm),Psamp(3*natm))
allocate(cov_R(3*natm,3*natm),cov_P(3*natm,3*natm),mean_R(3*natm),mean_P(3*natm),vec(3*natm),mat(3*natm,3*natm))
cov_R=0.0D0
cov_P=0.0D0
mean_R=0.0D0
mean_P=0.0D0



do i=1,nssmp
	
	if(i_typesamp == 0)call sample_classical(Rsamp,Psamp,eigvalarr,eigvecmat,atmmass,beta)
	if(i_typesamp == 1)call sample_wigner(Rsamp,Psamp,eigvalarr,eigvecmat,atmmass,beta)
	Rsamp=Rsamp+coord 

	if(i_typeout == 1)then 
		write(c10, '(I10)') i 
		call system("mkdir -p sample_"//trim(adjustl(c10)))
		open(11,file="sample_"//trim(adjustl(c10))//"/R0.xyz",status="replace")
		open(12,file="sample_"//trim(adjustl(c10))//"/P0.xyz",status="replace")
	end if

	write(11,*) natm
	write(11,*) "Sampled coordinates at step ",i
	do j=1,natm
		write(11,"(a2,3E18.8)") atom_names(int(atomid(j))), Rsamp((j-1)*3+1) * 0.529177249, Rsamp((j-1)*3+2)* 0.529177249, Rsamp((j-1)*3+3)* 0.529177249
	end do

	write(12,*) natm
	write(12,*) "Sampled coordinates at step ",i
	do j=1,natm
		write(12,"(a2,3E18.8)") atom_names(int(atomid(j))), Psamp((j-1)*3+1), Psamp((j-1)*3+2), Psamp((j-1)*3+3)
	end do

	if(i_typeout == 1)then 
		close(11)
		close(12)
	end if

	mean_R=mean_R+Rsamp
	mean_P=mean_P+Psamp
	call outer_product(Rsamp,Rsamp,mat)
	cov_R=cov_R+mat
	call outer_product(Psamp,Psamp,mat)
	cov_P=cov_P+mat

end do

if(i_typeout == 0)then 
	close(11)
	close(12)
end if

mean_R=mean_R/nssmp
mean_P=mean_P/nssmp
call outer_product(mean_R,mean_R,mat)
cov_R=cov_R/nssmp-mat 
call outer_product(mean_P,mean_P,mat)
cov_P=cov_P/nssmp-mat







! ====================================================================================

! call diagsymat(cov_R,eigvecmat,eigvalarr,istat)
! mat=0.0D0
! do i=1,3*natm
! 	mat(i,i)=1.0/eigvalarr(i)
! end do
! cov_R=matmul(matmul(eigvecmat,mat),transpose(eigvecmat))*beta 

! call diagsymat(cov_P,eigvecmat,eigvalarr,istat)
! mat=0.0D0
! do i=1,3*natm
! 	mat(i,i)=1.0/eigvalarr(i)
! end do
! cov_P=matmul(matmul(eigvecmat,mat),transpose(eigvecmat))*beta

! write(*,*) "For checking Sampling:"
! write(*,*) "Mean of sampled coordinates:"
! write(*,"(6(E18.8))") mean_R * 0.529177249 !Convert to Angstrom
! write(*,*) "Mean of sampled momentums:"
! write(*,"(6(E18.8))") mean_P !Momentum in a.u.
! write(*,*) "Covariance of sampled coordinates:"
! mat=matmul(Hess*beta,cov_R) 
! call showmatgau(mat,"Covariance of sampled coordinates")
! mat=0.0D0
! do i=1,3*natm
! 	mat(i,i)=1.0/massmat(i,i)
! end do
! mat=matmul(mat*beta,cov_P)
! write(*,*) "Covariance of sampled momentums:"
! call showmatgau(mat,"Covariance of sampled momentums")


! !Convert force constant to harmonic frequencies
! amu2kg=1.66053878D-27
! b2m=0.529177249D-10
! au2J=4.35974434D-18
! eigvalarr=eigvalarr*au2J/b2m**2/amu2kg !First convert force constant from a.u. to SI
! allocate(freq(nmodeall),wavenum(nmodeall))
! do i=1,nmodeall
! 	if (eigvalarr(i)<0) then
! 		freq(i)=-dsqrt(abs(eigvalarr(i)))/(2*pi)
! 	else
! 		freq(i)=dsqrt(eigvalarr(i))/(2*pi)
! 	end if
! end do

! !Sort all modes according to absolute value of frequencies (from low to high)
! allocate(tmpvec(nmodeall))
! do i=1,nmodeall
! 	do j=i+1,nmodeall
! 		if (abs(freq(i))>abs(freq(j))) then
! 			temp=freq(i)
! 			freq(i)=freq(j)
! 			freq(j)=temp
! 			tmpvec=eigvecmat(:,i)
! 			eigvecmat(:,i)=eigvecmat(:,j)
! 			eigvecmat(:,j)=tmpvec
! 		end if
! 	end do
! end do
! !Now the first six modes are considered as overall motions, we now sort remaining modes from low to high
! do i=7,nmodeall
! 	do j=i+1,nmodeall
! 		if (freq(i)>freq(j)) then
! 			temp=freq(i)
! 			freq(i)=freq(j)
! 			freq(j)=temp
! 			tmpvec=eigvecmat(:,i)
! 			eigvecmat(:,i)=eigvecmat(:,j)
! 			eigvecmat(:,j)=tmpvec
! 		end if
! 	end do
! end do
! wavenum=freq/2.99792458D10

! !Convert normal coordinates from mass-weighted basis to Cartesian basis
! allocate(normvec(nmodeall,nmodeall))
! do i=1,natm
! 	do j=(i-1)*3+1,i*3
! 		normvec(j,:)=eigvecmat(j,:)/dsqrt(atmmass(i))
! 	end do
! end do
! do i=1,nmodeall !Normalization
! 	tmp=dsqrt(sum(normvec(:,i)**2))
! 	normvec(:,i)=normvec(:,i)/tmp
! end do
! write(*,*)
! call showmatgau(normvec(:,7:),"Normal coordinates (columns)",form="1x,f9.4,4x")

! !Output final frequencies
! write(*,*)
! write(*,*) "The frequencies (cm-1) corresponding to overall translation and rotation:"
! write(*,"(6(f12.5))") wavenum(1:6)
! write(*,*) "Harmonic vibrational frequencies:"
! idx=0
! do i=7,nmodeall
! 	idx=idx+1
! 	write(*,"(' Mode',i5,':',E15.5,' Hz',f15.5,' cm-1')") idx,freq(i),wavenum(i)
! end do
! write(*,*)

if (isilent==0) then
	write(*,*) "Press ENTER to exit"
	read(*,*)
end if

contains

!!------------ Diagonalize a symmetry matrix 
!Repack the extremely complex "DSYEV" routine in lapack to terse form
!if istat/=0, means error occurs
subroutine diagsymat(mat,eigvecmat,eigvalarr,istat)
integer istat
real*8 mat(:,:),eigvecmat(:,:),eigvalarr(:)
real*8,allocatable :: lworkvec(:)
isize=size(mat,1)
allocate(lworkvec(3*isize-1))
call DSYEV('V','U',isize,mat,isize,eigvalarr,lworkvec,3*isize-1,istat)
eigvecmat=mat
mat=0D0
forall (i=1:isize) mat(i,i)=eigvalarr(i)
end subroutine

!!!-------- Find the line where the label first appears in fileid
!Return ifound=1 if found the label, else return 0
!Default is rewind, if irewind=0 will not rewind
!If current line already has the label, calling this subroutine will do nothing
subroutine loclabel(fileid,label,ifound,irewind)
integer fileid,ierror
integer,optional :: ifound,irewind
character*200 c200
CHARACTER(LEN=*) label
if ((.not.present(irewind)).or.(present(irewind).and.irewind==1)) rewind(fileid)
do while(.true.)
	read(fileid,"(a)",iostat=ierror) c200
	if (index(c200,label)/=0) then
		backspace(fileid)
		if (present(ifound)) ifound=1 !Found result
		return
	end if
	if (ierror/=0) exit
end do
if (present(ifound)) ifound=0
end subroutine

!!------ Display matrix similar to gaussian program, automatically switch to next screen
subroutine showmatgau(mat,label,insemi,form,fileid,useri1,useri2,inncol,titlechar)
!Number of columns is always 5, unadjustable
!"Label" is the title, if content is empty, title will not be printed
!If semi==1, only lower and diagonal element will be shown
!"form" is the format to show data, default is D14.6, can pass into such as "f14.8", total width should be 14 characters
!fildid is destination, 6 corresponds output to screen
!useri1 and useri2 is the dimension of the matrix, default or =-1 is determine automatically
!inncol seems controls spacing between number labels of each frame
!titlechar default is "i8,6x", if you manually set inncol, you also set this to broaden or narrow
implicit real*8(a-h,o-z)
real*8 :: mat(:,:)
character(*),optional :: label,form,titlechar
integer,optional :: insemi,fileid,useri1,useri2,inncol
integer :: semi,ides,ncol
semi=0
ides=6
ncol=5
i1=size(mat,1)
i2=size(mat,2)
if (present(useri1).and.useri1/=-1) i1=useri1
if (present(useri2).and.useri1/=-1) i2=useri2
if (present(insemi)) semi=insemi
if (present(fileid)) ides=fileid
if (present(inncol)) ncol=inncol
if (present(label).and.label/='') write(ides,*) "************ ",label," ************"
nt=ceiling(i2/float(ncol))
do i=1,nt !How many frames
	ns=(i-1)*5+1 !This frame starts from where
	if (i/=nt) ne=(i-1)*ncol+ncol !This frame end to where
	if (i==nt) ne=i2
	!Write basis number in separate line
	write(ides,"(6x)",advance='no')
	do j=ns,ne
		if (present(titlechar)) then
			write(ides,'('//titlechar//')',advance='no') j
		else
			write(ides,"(i8,6x)",advance='no') j
		end if
	end do
	write(ides,*)
	!Write content in each regular line
	do k=1,i1
		if (k<ns.and.semi==1) cycle !The lines have been outputted are skipped
		write(ides,"(i6)",advance='no') k
		do j=ns,ne
			if (semi==1.and.k<j) cycle !Upper trigonal element were passed
			if (present(form)) then
				write(ides,'('//form//')',advance='no') mat(k,j)			
			else
				write(ides,"(D14.6)",advance='no') mat(k,j)
			end if
		end do
		write(ides,*) !Change to next line
	end do
end do
end subroutine


subroutine box_muller(x1,x2,sigma,miu)
    implicit none
    !real*8, parameter :: pi = 3.14159265358979323846
    !real*8, parameter :: hbar    = 1.0
    real*8 x1,x2,u1,u2,sigma,miu
    call RANDOM_NUMBER(u1)
    call RANDOM_NUMBER(u2)
    x1=sqrt(-2*log(u1))*cos(2*pi*u2)
    x2=sqrt(-2*log(u1))*sin(2*pi*u2)
    !x1=(x1-miu)/sigma
    !x2=(x2-miu)/sigma
    x1=x1*sigma+miu
    x2=x2*sigma+miu
end subroutine



subroutine sample_classical(Rsamp,Psamp,freq2,U,mass,beta)
implicit none
real*8 Rsamp(:),Psamp(:),freq2(:),U(:,:),mass(:),beta
real*8,allocatable :: mass3n(:)
integer :: ire,nsize
real*8 :: x2, cutoff,au2wn 


au2wn = 219474.6313702
cutoff=100 / au2wn

nsize=size(Rsamp,1)
allocate(mass3n(nsize))
do ire=1,nsize/3
	mass3n((ire-1)*3+1)=mass(ire)
	mass3n((ire-1)*3+2)=mass(ire)
	mass3n((ire-1)*3+3)=mass(ire)
end do

do ire=1,nsize
	!Sample a random number from normal distribution
	if (freq2(ire)<=cutoff**2) then
		Rsamp(ire)=0.0D0 !If frequency is zero or negative, set position to zero
		Psamp(ire)=0.0D0 !Set momentum to zero
	else
		call box_muller(Rsamp(ire),x2,sqrt(1/(beta*freq2(ire))),0.0D0)
		call box_muller(Psamp(ire),x2,sqrt(1/(beta)),0.0D0)

	end if
	
end do

Rsamp=matmul(U,Rsamp)/dsqrt(mass3n) !Convert to Cartesian coordinates
Psamp=matmul(U,Psamp)*dsqrt(mass3n) !Convert to Cartesian coordinates

deallocate(mass3n)
end subroutine



subroutine sample_wigner(Rsamp,Psamp,freq2,U,mass,beta)
implicit none
real*8 Rsamp(:),Psamp(:),freq2(:),U(:,:),mass(:),beta
real*8,allocatable :: mass3n(:),Qcor(:)
integer :: ire,nsize
real*8 :: x2, cutoff,cutoff2
real*8 :: amu2kg,b2m,au2J,au2wn 



au2wn = 219474.6313702

cutoff=0 / au2wn 
cutoff2=30  / au2wn

nsize=size(Rsamp,1)
allocate(mass3n(nsize))
do ire=1,nsize/3
	mass3n((ire-1)*3+1)=mass(ire)
	mass3n((ire-1)*3+2)=mass(ire)
	mass3n((ire-1)*3+3)=mass(ire)
end do


allocate(Qcor(nsize))
Qcor=0.0D0
do ire=1,nsize
	if (freq2(ire)<=0) then
		Qcor(ire)=tanh(0.5*beta*sqrt(-freq2(ire))) /(0.5*beta*sqrt(-freq2(ire)))
	else
		Qcor(ire)=0.5*beta*sqrt(freq2(ire))/tanh(0.5*beta*sqrt(freq2(ire))) 
	end if
	
end do

do ire=1,nsize
	!Sample a random number from normal distribution
	! print*,freq2(ire),cutoff**2
	if (freq2(ire)<=cutoff**2) then
		
		Rsamp(ire)=0.0D0 !If frequency is zero or negative, set position to zero
		Psamp(ire)=0.0D0 !Set momentum to zero
	else if (freq2(ire)>cutoff**2 .and. freq2(ire)<=cutoff2**2) then
		
		call box_muller(Rsamp(ire),x2,sqrt(1/(beta*freq2(ire))),0.0D0)
		call box_muller(Psamp(ire),x2,sqrt(1/(beta)),0.0D0)
	else
		
		if(beta>1e5)then
			call box_muller(Rsamp(ire),x2,sqrt(0.5/(sqrt(freq2(ire)))),0.0D0)
			call box_muller(Psamp(ire),x2,sqrt(0.5 * sqrt(freq2(ire))),0.0D0)
		else
			call box_muller(Rsamp(ire),x2,sqrt(Qcor(ire)/(beta*freq2(ire))),0.0D0)
			call box_muller(Psamp(ire),x2,sqrt(Qcor(ire)/(beta)),0.0D0)
		end if
	end if
	
end do
! stop
Rsamp=matmul(U,Rsamp)/dsqrt(mass3n) !Convert to Cartesian coordinates
Psamp=matmul(U,Psamp)*dsqrt(mass3n) !Convert to Cartesian coordinates

deallocate(mass3n)
end subroutine

subroutine outer_product(a,b,c)
implicit none
real*8, intent(in) :: a(:), b(:)
real*8, intent(out) :: c(:,:)
integer :: i, j, n, m
n = size(a)
m = size(b)
c = 0.0D0
do i = 1, n
	do j = 1, m
		c(i,j) = a(i) * b(j)
	end do
end do

end subroutine 


end program
