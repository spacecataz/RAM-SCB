!============================================================================
!    Copyright (c) 2016, Los Alamos National Laboratory
!    All rights reserved.
!============================================================================

MODULE ModRamInit
! Contains subroutines for initialization of RAM

  implicit none

  contains
!==============================================================================
  subroutine ram_allocate
    ! Allocate arrays for variables used in RAM calculations

    use ModRamVariables ! Need to allocate and initialize all the variables
    use ModRamParams,    ONLY: IsComponent, FixedComposition
    use ModRamSpecies,   ONLY: RAMSpecies, nSpecies
    use ModRamGrids,     ONLY: nR, nRExtend, nT, nE, nPa, &
                               Slen, ENG, NCF, NL, nS, nX, &
                               RadMaxScb, RadiusMin, RadiusMax, &
                               NameVar, NCF_emic, ENG_emic
    use ModScbGrids,     ONLY: nthe

    use ModRamIO,        ONLY: read_CHEX_file

    implicit none

    character(len=100) :: tempStr
    character(len=2)   :: tempSpc
    integer :: i, iS, nChar

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='ram_allocate'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Determine max size of grid for full SCB coverage via SWMF:
    nRExtend = (nR-1)*(RadMaxScb-RadiusMin)/(RadiusMax-RadiusMin) + 1
    nX = NPA

    if(DoTest)write(*,*)'IM/',NameSub,': nRextend = ', nRextend

    ! Initialize species being run
    ! If needed, read the charge exchange cross sections file
    ! If using a fixed composition model, set the compositions as well
    allocate(species(nS))
    tempStr = trim(NameVar)
    do i = 1, nS
       nChar = index(tempStr, ' ')
       if (nChar.eq.-1) nChar = len(trim(tempStr))
       tempSpc = trim(tempStr(1:nChar))
       do iS = 1, nSpecies
          if (tempSpc == RAMSpecies(iS)%s_code) then
             species(i) = RAMSpecies(iS)
             if (trim(RAMSpecies(iS)%CEX_file).ne.'na') then
                call read_CHEX_file(trim(RAMSpecies(iS)%CEX_file), species(i))
             endif
             if (FixedComposition) species(i)%s_comp = composition(i)
          endif
       enddo
       tempStr = trim(tempStr(nChar+1:len(tempStr)))
    enddo

  !!!!!!!! Allocate Arrays
  ! Initialization of variables is not strictly needed for modern compilers,
  ! but it removes any potential ambiguity
    ALLOCATE(outsideMGNP(nR,nT))
    outsideMGNP = 0

  ! Main RAM Variables
    ALLOCATE(F2(NS,NR,NT,NE,NPA), FLUX(NS,NR,NT,NE,NPA),PPerH(NR,NT), PParH(NR,NT), &
             PPerE(NR,NT), PParE(NR,NT), PPerO(NR,NT),PParO(NR,NT), PPerHe(NR,NT), &
             PParHe(NR,NT), PAllSum(NR,NT), PParSum(NR,NT), PPerT(NS,NR,NT), &
             PParT(NS,NR,NT), FNHS(NR+1,NT,NPA), FNIS(NR+1,NT,NPA), BOUNHS(NR+1,NT,NPA), &
             BOUNIS(NR+1,NT,NPA), dIdt(NR+1,NT,NPA), dBdt(NR+1,NT), dIbndt(NR+1,NT,NPA), &
             HDNS(NR+1,NT,NPA), BNES(NR+1,NT), dHdt(nR+1,nT,nPa), ODNS(NR+1,NT,NPA), NDNS(NR+1,NT,NPA))
    F2 = 0._dp; FLUX = 0._dp; PPerH = 0._dp; PParH = 0._dp; PPerE = 0._dp; PParE = 0._dp; PPerO = 0._dp
    PParO = 0._dp; PPerHe = 0._dp; PParHe = 0._dp; PAllSum = 0._dp; PParSum = 0._dp; PPerT = 0._dp
    PParT = 0._dp; FNHS = 0._dp; FNIS = 0._dp; BOUNHS = 0._dp; BOUNIS = 0._dp; dIdt = 0._dp
    dBdt = 0._dp; dIbndt = 0._dp; HDNS = 0._dp; BNES = 0._dp; dHdt = 0._dp
    ODNS = 0._dp; NDNS = 0._dp

  ! ModRamPlasmasphere Variables
    ALLOCATE(flux_volume(nR,nT), tau(nR,nT), uL(nR,nT), uT(nT,nR), smLon(nR,nT), smLat(nR,nT), &
             nECR(nR,nT), xRAM(nthe,nR,nT), yRAM(nthe,nR,nT), zRAM(nthe,nR,nT))

  ! ModRamCouple and IM_wrapper Variables
    ALLOCATE(NAllSum(nR,nT), DensO(nR,nT), DensH(nR,nT), DensHe(nR,nT), HPAllSum(nR,nT), &
             OPAllSum(nR,nT), HePAllSum(nR,nT), ePAllSum(nR,nT), HNAllSum(nR,nT), &
             ONAllSum(nR,nT), HeNAllSum(nR,nT))
    NAllSum = -1._dp; DensO = -1._dp; DensH = -1._dp; DensHe = -1._dp; HPAllSum = -1._dp
    OPAllSum = -1._dp; HePAllSum = -1._dp; ePAllSum = -1._dp; HNAllSum = -1._dp
    ONAllSum = -1._dp; HeNAllSum = -1._dp

  ! ModRamInit Variables
    ALLOCATE(RMAS(NS), V(NS,NE), VBND(NS,NE), GREL(NS,NE), GRBND(NS,NE), FACGR(NS,NE), &
             EPP(NS,NE), ERNH(NS,NE), UPA(NR), WE(NE), DE(NE), EKEV(NE), EBND(NE), &
             PHI(NT), LT(NT), MLT(NT), MU(NPA), DMU(NPA), WMU(NPA), PAbn(NPA), LZ(NR+1), &
             RLZ(NR+1), AMLA(Slen), BE(NR+1,Slen), GridExtend(NRExtend), ZRPabn(NR,NPA,Slen), &
             FFACTOR(NS,NR,NE,NPA), PA(NPA))
    RMAS = 0._dp; V = 0._dp; VBND = 0._dp; GREL = 0._dp; GRBND = 0._dp; FACGR = 0._dp; EPP = 0._dp
    ERNH = 0._dp; UPA = 0._dp; WE = 0._dp; DE = 0._dp; EKEV = 0._dp; EBND = 0._dp; PHI = 0._dp
    LT = 0._dp; MLT = 0._dp; MU = 0._dp; DMU = 0._dp; WMU = 0._dp; PAbn = 0._dp; LZ = 0._dp; PA = 0._dp
    RLZ = 0._dp; AMLA = 0._dp; BE = 0._dp; GridExtend = 0._dp; ZrPabn = 0._dp; FFACTOR = 0._dp

  ! ModRamWPI Variables
    ALLOCATE(WALOS1(NR,NE), WALOS2(NR,NE), WALOS3(NR,NE), fpofc(NCF), NDVVJ(NR,ENG,NPA,NCF), &
             NDAAJ(NR,ENG,NPA,NCF), ENOR(ENG), ECHOR(ENG), BDAAR(NR,NT,ENG,NPA), &
             CDAAR(NR,NT,NE,NPA))
    WALOS1 = 0._dp; WALOS2 = 0._dp; WALOS3 = 0._dp; fpofc = 0._dp; NDVVJ = 0._dp; NDAAJ = 0._dp
    ENOR = 0._dp; ECHOR = 0._dp; BDAAR = 0._dp; CDAAR = 0._dp

    ALLOCATE(Daa_emic_h(NR,ENG_emic, NPA, NCF_emic), Daa_emic_he(NR,ENG_emic, NPA, NCF_emic), &
             EKEV_emic(ENG_emic), fp2c_emic(NCF_emic), Ihs_emic(4,20,25), Ihes_emic(4,20,25))
    Daa_emic_h = 0.0; Daa_emic_he = 0.0; EKEV_emic = 0.0; fp2c_emic = 0.0
    Ihs_emic = 0.0; Ihes_emic = 0.0

  ! ModRamLoss Variables
    ALLOCATE(ATLOS(nS,NR,NE), CHARGE(nS,NR,NT,NE,NPA), COULE(nS,nE,nPa), COULI(nS,nE,nPa), &
             ATA(nS,nE,nPa), GTA(nS,nE,nPa), GTAE(nS,nE,nPa), GTAI(nS,nE,nPa), &
             CEDR(nS,nE,nPa), CIDR(nS,nE,nPa), FLC_coef(4,NR,NT,NE,NPA), &
             r_curvEq(NR,NT), zeta1Eq(NR,NT), zeta2Eq(NR,nT))
    ATLOS = 0._dp; CHARGE = 0._dp; COULE = 0.0; COULI = 0.0; ATA = 0.0; GTA = 0.0
    GTAE = 0.0; GTAI = 0.0; CEDR = 0.0; CIDR = 0.0; FLC_coef = 0._dp
    r_curvEq = 0.0; zeta1Eq = 0.0; zeta2Eq = 0.0

  ! ModRamEField Variables
    ALLOCATE(VT(NR+1,NT), EIR(NR+1,NT), EIP(NR+1,NT), VTOL(NR+1,NT), VTN(NR+1,NT))
    VT = 0._dp; EIR = 0._dp; EIP = 0._dp; VTOL = 0._dp; VTN = 0._dp

  ! ModRamBoundary Variables
    ALLOCATE(FGEOS(NS,NT,NE,NPA))
    FGEOS = 0._dp

  ! ModRamDrift Variables
    ALLOCATE(DtDriftR(nS), DtDriftP(nS), DtDriftE(nS), DtDriftMu(nS))
    DtDriftR = 0._dp; DtDriftP = 0._dp; DtDriftE = 0._dp; DtDriftMu = 0._dp

  ! ModRamRun Variables
    ALLOCATE(SETRC(NS), ELORC(NS), LSDR(NS), LSCHA(NS), LSATM(NS), LSCOE(NS), &
             LSCSC(NS), LSWAE(NS), XNN(NS,NR), XND(NS,NR), LNCN(NS,NR), LNCD(NS,NR), &
             LECN(NS,NR), LECD(NS,NR), ENERN(NS,NR), ENERD(NS,NR), ATEW(NR,NT,NE,NPA), &
             ATAW(NR,NT,NE,NPA), ATAC(NR,NT,NE,NPA), ATEC(NR,NT,NE,NPA), XNE(NR,NT), &
             ATMC(NR,NT,NE,NPA), ATAW_emic_h(NR,NT,NE,NPA), ATAW_emic_he(NR,NT,NE,NPA),&
             ESUM(NS), NSUM(NS))
    SETRC = 0._dp; ELORC = 0._dp; LSDR = 0._dp; LSCHA = 0._dp; LSATM = 0._dp; LSCOE = 0._dp
    LSCSC = 0._dp; LSWAE = 0._dp; XNN = 0._dp; XND = 0._dp; LNCN = 0._dp; LNCD = 0._dp
    LECN = 0._dp; LECD = 0._dp; ENERN = 0._dp; ENERD = 0._dp; ATEW = 0._dp; ATAW = 0._dp
    ATAC = 0._dp; ATEC = 0._dp; XNE = 0._dp; ATMC = 0._dp; ATAW_emic_h = 0._dp; ATAW_emic_he = 0._dp;
    NECR = 0._dp; ESUM = 0._dp; NSUM = 0._dp
  !!!!!!!!!

  end subroutine ram_allocate

!==================================================================================================
  subroutine ram_deallocate
    ! Deallocate all allocated arrays for cleanup

    use ModRamVariables ! Need to deallocate all variables

    implicit none

  !!!!!!!! Deallocate Arrays
    DEALLOCATE(outsideMGNP)
    DEALLOCATE(species)
  ! Main RAM Variables
    DEALLOCATE(F2, FLUX,PPerH, PParH, PPerE, PParE, PPerO, PParO, PPerHe, PParHe, &
               PAllSum, PParSum, PPerT, PParT, FNHS, FNIS, BOUNHS, BOUNIS, dIdt, &
               dBdt, dIbndt, HDNS, BNES, dHdt)
  ! ModRamPlasmasphere Variables
    DEALLOCATE(flux_volume, tau, uL, uT, smLon, smLat, NECR, xRAM, yRAM, zRAM)
  ! ModRamInit Variables
    DEALLOCATE(RMAS, V, VBND, GREL, GRBND, FACGR, EPP, ERNH, UPA, WE, DE, EKEV, &
               EBND, PHI, LT, MLT, MU, DMU, WMU, PAbn, LZ, RLZ, AMLA, BE, GridExtend, &
               ZRPabn, FFACTOR)
  ! ModRamWPI Variables
    DEALLOCATE(WALOS1, WALOS2, WALOS3, fpofc, NDVVJ, NDAAJ, ENOR, ECHOR, BDAAR, &
               CDAAR)
    DEALLOCATE(Daa_emic_h, Daa_emic_he, EKEV_emic, fp2c_emic, Ihs_emic, Ihes_emic)
  ! ModRamLoss Variables
  !  DEALLOCATE(ATLOS, ACHAR)
    DEALLOCATE(r_curvEq, zeta1Eq, zeta2Eq, FLC_coef)
  ! ModRamEField Variables
    DEALLOCATE(VT, EIR, EIP, VTOL, VTN)
  ! ModRamBoundary Variables
    DEALLOCATE(FGEOS)
  ! ModRamDrift Variables
    DEALLOCATE(DtDriftR, DtDriftP, DtDriftE, DtDriftMu)
  ! ModRamRun Variables
    DEALLOCATE(SETRC, ELORC, LSDR, LSCHA, LSATM, LSCOE, LSCSC, LSWAE, XNN, XND, &
               LNCN, LNCD, LECN, LECD, ENERN, ENERD, ATEW, ATAW, ATAC, ATEC, &
               XNE, ATMC, ATAW_emic_h, ATAW_emic_he,ESUM, NSUM)
  !!!!!!!!!


  end subroutine ram_deallocate

!==================================================================================================
  SUBROUTINE ram_init
    ! Initialize RAM

    !!!! Module Variables
    use ModRamParams,    ONLY: DoUseWPI, DoUseBASDiff, IsRestart, IsComponent,DoUseEMIC
    use ModRamMain,      ONLY: DP, PathRestartIn, nIter
    use ModRamTiming,    ONLY: TimeRamStart, TimeMax, TimeRamRealStart, TimeRamNow, &
                               TimeRamElapsed, TimeMax, TimeRestart, TimeRamFinish, &
                               TOld
    use ModRamGrids,     ONLY: RadiusMax, RadiusMin, nR, nRExtend, nT, dR, dPhi, nS
    use ModRamVariables, ONLY: PParH, PPerH, PParHe, PPerHe, PParO, PPerO, PParE, &
                               PPerE, LSDR, LSCHA, LSATM, LSCOE, LSCSC, LSWAE, ELORC, &
                               SETRC, XNN, XND, ENERN, ENERD, LNCN, LNCD, LECN, LECD, &
                               Lz, GridExtend, Phi, kp, F107, AE, species, ESUM, NSUM
    use ModScbVariables, ONLY: radRaw, azimRaw
    !!!! Modules Subroutines/Functions
    use ModRamWPI,     ONLY: WAPARA_HISS, WAPARA_BAS, WAPARA_CHORUS, WAVEPARA1, WAVEPARA2,&
                             WAPARA_EMIC
    use ModRamIndices, ONLY: init_indices, get_indices
    !!!! Share Modules
    use ModTimeConvert, ONLY: TimeType, time_real_to_int, time_int_to_real
    use ModNumConst,    ONLY: cTwoPi
    use ModIOUnit,      ONLY: UNITTMP_

    implicit none

    type(timetype) :: TimeRamStop

    real(DP) :: dPh

    integer :: iR, iPhi, j, k, iS
    integer :: nrIn, ntIn, neIn, npaIn
    logical :: TempLogical
    logical :: StopCommand, IsStopTimeSet
    character(len=100) :: StringLine, NameCommand, RestartFile

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='ram_init'
    !------------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! Set up time information for the simulation
    if (IsRestart) then
       RestartFile=PathRestartIn//'/restart_info.txt'
       open(unit=UnitTMP_, file=trim(RestartFile), status='old')
       read(UnitTMP_,*)StringLine
       read(UnitTMP_, '(a25,i4.4, 2i2.2, 1x, 3i2.2)')StringLine, &
            TimeRamStart%iYear, TimeRamStart%iMonth, TimeRamStart%iDay, &
            TimeRamStart%iHour, TimeRamStart%iMinute, TimeRamStart%iSecond
       TimeRamStart%FracSecond=0.0
       read(UnitTMP_,'(a25, f15.4)') StringLine, TimeRestart
       read(UnitTMP_,'(a25, i15)') StringLine, nIter
       read(UnitTMP_, *) StringLine
       read(UnitTMP_, '(a25, 4i3)') StringLine, nrIn, ntIn, neIn, npaIn
       close(UnitTMP_)
       call time_int_to_real(TimeRamStart)
       TimeRamRealStart%Time = TimeRamStart%Time + TimeRestart
       TimeRamElapsed = TimeRestart
       call time_real_to_int(TimeRamRealStart)
    else
       TimeRamElapsed = 0
       TimeRamRealStart = TimeRamStart
       TimeRestart = 0
    end if
    TimeRamNow = TimeRamRealStart
    TOld = TimeRamElapsed

    ! Calculate TimeMax
    if (IsComponent) then
       TimeRamNow = TimeRamRealStart
    else
       If (abs(TimeMax).le.1e-9) then
          TimeMax = TimeRamFinish%Time-TimeRamStart%Time
       else
          TimeMax = TimeRamElapsed + TimeMax
       endif
    endif

    TimeRamStop%Time = TimeRamStart%Time + TimeMax
    call time_real_to_int(TimeRamStop)
    call init_indices(TimeRamRealStart, TimeRamStop)
    call get_indices(TimeRamNow%Time, Kp, f107, AE)

  !!!!!!!!! Zero Values
    ! Initial loss is zero
    LNCN  = 0._dp; LNCD  = 0._dp
    LECN  = 0._dp; LECD  = 0._dp
    LSDR  = 0._dp
    LSCHA = 0._dp
    LSATM = 0._dp
    LSCOE = 0._dp
    LSCSC = 0._dp
    LSWAE = 0._dp
    ELORC = 0._dp
    SETRC = 0._dp

    ! Initial energy and density
    XNN   = 0._dp; XND   = 0._dp; NSUM = 0._dp
    ENERN = 0._dp; ENERD = 0._dp; ESUM = 0._dp
  !!!!!!!!!

  !!!!!!!!!! Initialize grid.
    ! Radial distance
    dR = (RadiusMax - RadiusMin)/(nR - 1)
    do iR = 1, nR+1 ! DANGER WE SHOULD CHANGE THIS AND ALL ARRAYS USING NR+1
       Lz(iR) = RadiusMin + (iR - 1)*dR
    end do

    ! Create extended radial grid for coupling:
    do iR=1, nRextend
       GridExtend(iR) = RadiusMin + (iR-1)*dR
    end do

    ! Longitude in radians
    dPh = cTwoPi/(nT - 1)
    do iPhi = 1, nT
       Phi(iPhi) = (iPhi - 1)*dPh
    end do

    if(DoTest)then
       write(*,*)'IM: GRID INITIALIZED'
       write(*,*)'  nT, nR, nRextend = ', nT, nR, nRextend
       write(*,*)'  dPh, dR = ', dPh, dR
       write(*,*)'  GridExtend = ', GridExtend
       write(*,*)'  Lz (Ram-only grid) = ', Lz
       write(*,*)'  Phi = ', phi
    end if

    ! Arrays that are needed in the SCB calculations.
    DO j = 0,nR
       radRaw(j) = RadiusMin + ((RadiusMax+dR)-RadiusMin) * REAL(j,DP)/REAL(nR,DP)
    END DO
    DO k = 1,nT
       azimRaw(k) = 24.0 * REAL(k-1,DP)/REAL(nT-1,DP)
    END DO

    ! Intialize arrays needed for RAM
    do iS=1,nS
       call Arrays(iS)
       ! If needed, initialize wave particle interaction arrays for the given species
       if (species(iS)%WPI) then
          IF (DoUseWPI) THEN
             CALL WAPARA_HISS(iS)
             IF (DoUseBASdiff) then
                CALL WAPARA_BAS(iS)
             ELSE
                CALL WAPARA_CHORUS(iS)
             ENDIF
          ELSE
             CALL WAVEPARA1(iS)
             CALL WAVEPARA2(iS)
          end if
       ENDIF
    end do
    if(DoUseEMIC) call WAPARA_EMIC ! call only once; the Daa is same for each species

  END SUBROUTINE ram_init

!**************************************************************************
!                               ARRAYS
!                       Set up all the arrays
!**************************************************************************
  SUBROUTINE ARRAYS(S)
    !!!! Module Variables
    use ModRamMain,      ONLY: DP
    use ModRamConst,     ONLY: RE, PI, MP, CS, Q, HMIN
    use ModRamGrids,     ONLY: RadiusMax, RadiusMin, NR, NPA, Slen, NT, NE, &
                               NLT, EnergyMin, dR, dPhi, nS
    use ModRamParams,    ONLY: DoUsePlasmasphere
    use ModRamVariables, ONLY: amla, DL1, Lz, RLz, IR1, EKEV, Mu, WMu, DMu, &
                               RMAS, WE, DE, EBND, GRBND, V, Pa, Pabn, UPA, &
                               FFACTOR, GREL, ZrPabn, VBND, PHI, BE, MLT, &
                               ERNH, EPP, FACGR, CONF1, CONF2, IP1, &
                               MDR, RFACTOR, species
    !!!! Module Subroutines/Functions
    use ModRamFunctions, ONLY: ACOSD, ASIND, COSD, SIND

    implicit none
    integer, intent(in) :: S

    real(DP) :: degrad, camlra, elb, rw, rwu
    real(DP) :: clc, spa, MUBOUN
    real(DP), ALLOCATABLE :: CONE(:),RLAMBDA(:)

    integer :: i, j, k, l, iml, ic, ip, iS

    ! Give Cone an array size of 4 more then the number of grid points in the
    ! radial direction.
    ! Give RLAMBDA an array size eaqual to the number of grid points in the
    ! picth angle dimension
    ALLOCATE(CONE(NR+4),RLAMBDA(NPA))
    CONE = 0.0; RLAMBDA = 0.0

    ! Grid size of L shell (in Earth Radii)
    DL1 = (RadiusMax - RadiusMin)/(nR - 1)

    degrad=pi/180. ! constant which converts degrees to radians
    amla(1)=0. ! Magnetic latitude grid in degrees
    DO I=2,6
      amla(i)=amla(i-1)+0.2
    ENDDO
    DO I=7,24
      amla(i)=amla(i-1)+0.5
    ENDDO
    DO I=25,Slen ! Slen is 55 defined in ModRamGrids.f90
      amla(i)=amla(i-1)+2.
    ENDDO
    ! amla produces an array which goes form 0 to 72 as of commit
    ! 10b4801d2d89e6fd663e248a46ecc4d65dd6fe11 (#368)
    ! on the main LANL fork.
    ! CBW note: Why does it stop short of 90 if it is suppose to be an array of
    ! magnetic latitude?
    !write(*,*) "amla:" ! debug print statements.
    !write(*,*) amla

    IR1=DL1/0.25                ! IR1 is 1/4 the size of the radial step
    MDR=DL1*RE                  ! Grid size (in meters) for Z=RO
    DO I=1,NR+1
      LZ(I)=RadiusMin+(I-1)*DL1 ! Grid L-Shells, there is a note in ram_init
      ! that LZ defined simiarlly there (NR+1) is dangerous and should be changed.
      RLZ(I)=RE*LZ(I)           ! Grid radial distance (from center of earth)
      DO IML=1,Slen
        camlra=amla(iml)*degrad ! radians of magnetic latitude.
        BE(I,IML)=0.32/LZ(I)**3*SQRT(1.+3.*SIN(camlra)**2)/COS(camlra)**6
      ENDDO
    END DO
    !write(*,*) "LZ:"
    !write(*,*) LZ
    !write(*,*) "RLZ:"
    !write(*,*) RLZ
    ! BE is the magnetic feild of earth in the dipole approximation depending on
    ! L-shell and equatorial pitch angle.

    DPHI=2.*PI/(NT-1)      ! Grid size for local time [rad]
    DO J=1,NT
      PHI(J)=(J-1)*DPHI    ! Magnetic local time in radian
      MLT(J)=PHI(J)*12./PI ! Magnetic local time in hour
    END DO
    IP1=(MLT(2)-MLT(1))/0.5 ! CBW note: IP1 is called in by two different
    ! modules, but is used only in this line. Must be defined here for use elsewhere

    DO iS=1,nS
      RMAS(iS)=MP*species(iS)%s_mass ! rest mass of each species (kg)
    END DO

    ! Calculate Kinetic Energy EKEV [keV] at center of energy bin, RW depends on
    ! NE. Set the values of WE and RW based on initial value of ELB
    ELB=EnergyMin ! Lower limit of energy in keV
    IF (abs(ELB-0.01).le.1e-9) THEN
      WE(1)=2.8E-3 !  |_._|___.___|____.____|______.______|
      RW=1.36      !    .     <   DE   >    <      WE     >
    END IF         !   EKEV                EBND
    IF (abs(ELB-0.1).le.1e-9) THEN ! relativistic
      WE(1)=3E-2
      RW=1.27
    END IF
    IF (abs(ELB-1.0).le.1e-9) THEN
      WE(1)=0.31
      RW=1.16
    END IF

    EKEV(1)=ELB+0.5*WE(1) ! set start of energy array in KeV based on
    ! scale determined above. Q, Elementary Charge in Columb, CS is the speed of light,
    ! RMAS array of size species-in-simulation contains rest mass for each species.
    GREL(S,1)=1.+EKEV(1)*1000.*Q/RMAS(S)/CS/CS ! GREL has units of kC
    V(S,1)=CS*SQRT(GREL(S,1)**2-1.)/GREL(S,1) ! velocity
    EBND(1)=ELB+WE(1) !
    GRBND(S,1)=1.+EBND(1)*1000.*Q/RMAS(S)/CS/CS
    VBND(S,1)=CS*SQRT(GRBND(S,1)**2-1.)/GRBND(S,1)
    DO K=1,NE-1
      WE(K+1)=WE(K)*RW            ! WE(K) [keV] is a power series
      EBND(K+1)=EBND(K)+WE(K+1)   ! E[keV] at bound of grid
      DE(K)=0.5*(WE(K)+WE(K+1))   ! Average of WE(K) and WE(K+1)
      EKEV(K+1)=EKEV(K)+DE(K)     ! E[keV] at center of grid
      GREL(S,K+1)=1.+EKEV(K+1)*1000.*Q/RMAS(S)/CS/CS
      V(S,K+1)=CS*SQRT(GREL(S,K+1)**2-1.)/GREL(S,K+1)      ! Veloc [m/s] at cent
      GRBND(S,K+1)=1.+EBND(K+1)*1000.*Q/RMAS(S)/CS/CS
      VBND(S,K+1)=CS*SQRT(GRBND(S,K+1)**2-1.)/GRBND(S,K+1) ! Veloc [m/s] at bound
    END DO
    DE(NE)=0.5*WE(NE)*(1.+RW) ! manually handel last entry of DE

    ! CONE - pitch angle loss cone in degree
    ! Dipole loss cone. This needs to be changed in future version to support a
    ! non-dipole loss cone -ME
    DO I=1,NR
      CLC=(RE+HMIN)/RLZ(I) ! HMIN is altitude of dense atmosphere (m), RLZ is
      ! Grid radial distance (m), RE is earth radius in (m). CLC is inverse L-Shell
      ! scaled by the raito of earth's radius to the radius of earth and its
      ! dense atmosphere.
      CONE(I)=ASIND(SQRT(CLC**3/SQRT(4.-3.*CLC))) ! Becuase of this angle values
      ! now depend on radius resolution!
    END DO
    !write(*,*) "Cone: Pre near boundary assignment"
    !write(*,*) CONE
    CONE(NR+1)=2.5 ! to calcul PA grid near 0 deg
    CONE(NR+2)=1.5
    CONE(NR+3)=1.
    CONE(NR+4)=0. ! THIS GERENERATES A JUMP IN CONE! with current LANL values
    ! CONE jumps back up from a lower value.
    !write(*,*) "Cone: Post near boundary assignment"
    !write(*,*) CONE

    ! PA is equatorial pitch angle in deg - PA(1)=90, PA(NPA)=0.
    ! recall: NPA is the grid points in pitch angle dimension
    ! MU is cosine of equatorial PA
                         ! |_._|___.___|____.____|______.______|
                         !   MU    <  DMU   >    <     WMU     >

    ! Goal of this if/else is to create several arrays: PA, and MU are explained
    ! above, RWU is a fixed number used for scaling the WMU array (found only in
    ! else portion), WMU same size as PA and MU, however contains explained
    ! below. DMU(L) is the amount that MU must change from MU(L) to MU(L+1)
    ! s.t. after NPA-1 steps MU has gone from 0 to 1.
    if (nPa == 90) then
      !write(*,*) "We are in If"
      ! Define the first and last entires of the arrays PA and MU manually.
      PA(1)=90.
      MU(1)=0.
      PA(NPA)=0.
      MU(NPA)=1.
      ! An unused scaling factor in the if case.
      RWU=0.98
      ! WMU(1) is the resolution in cos of the array MU
      ! i.e. WMU(L) is the value which should be added to MU(L) to produce MU(L+1)
      WMU(1)=(MU(NPA)-MU(1))/nPa ! the entire WMU array is filled with this value.
      ! This loop fills the remainder of the array WMU, DMU, MU, PA, and
      ! defines the boundary values of MU and PA in the vars MUBOUN and Pabn
      ! respectfully.
      DO L = 1, nPa-1
        WMU(L+1)=WMU(L) ! WMU is filled with a single value
        DMU(L)=0.5*(WMU(L)+WMU(L+1)) ! Fill DMU with needed step of MU
        MU(L+1)=MU(L)+DMU(L) ! apply needed step to MU(L) to reach MU(L+1)
        PA(L+1)=ACOSD(MU(L+1)) ! ACOSD(MU(L+1)) : PA :: DMU : MU
        MUBOUN = Mu(L) + 0.5*WMu(L) ! The value of MU on the boundary is
        ! the value of MU next to the boundary + .5 the resolution of steps in
        ! cos(angle).
        Pabn(L) = ACOSD(MUBOUN) ! Despite PA(NPA) being defined above to have
        ! PA(NPA) = 0, this DO loop will change PA(NPA) != 0. PAbn has a value
        ! between PA(L) and PA(L+1). Therefore PAbn(NPA) = 0 (as it takes the last
        ! half step that PA does not). Although PAbn(0) != 90 as PA(0) does.
      END DO
      DMU(nPa) = 0.5*(wMu(nPa)+wMu(nPa-1)) ! As the above loop cannot fill in
      ! DMU as WMU(NPA+1) does not exist, its final value is filled by hand
      ! using a reflection method for extending the array.

      Pabn(nPa) = 0. ! Although the above loop should ensure PAbn(NPA) = 0, it
      ! is defined to be so here.
    else
      !write(*,*) "We are in Else"
      ! explictly set the boundary of PA and MU
      PA(1)=90.
      MU(1)=0.
      PA(NPA)=0.
      MU(NPA)=1.
      ! RWU is a scaling factor applied to WMU so that the step from cos(angle)
      ! is slightly less then what it should be.
      RWU=0.98
      WMU(1)=(MU(NPA)-MU(1))/32 ! Array same size as PA
      !write(*,*) "PA is about to call Acosd for first portion of array"
      DO L=1,46
        WMU(L+1)=WMU(L)*RWU
        DMU(L)=0.5*(WMU(L)+WMU(L+1))
        MU(L+1)=MU(L)+DMU(L)
        PA(L+1)=ACOSD(MU(L+1))
      END DO
      !write(*,*) "PA called Acosd with no crash"
      ! Why is this being hardcoded, it effictively makes using an nPA less then
      ! 49 impossible! As it makes using NPA != 72 difficult as these numbers are
      ! adjusted for that value.
      PA(48)=18.7 ! CBW note: This imples that the 48th step from 90 to a step
      ! sort of 0 with 71 total steps is 18.7, but checking in python the 47th
      ! entry (python is 0 based indexing) in an arry similar to PA is 31.25.
      ! i.e. I think that assigning 18.7 to PA(48) is wrong.
      MU(48)=COSD(PA(48)) ! manually assign MU(48) based on value of PA(48)
      DMU(47)=(MU(48)-MU(47)) ! manually assign DMU(47) based on the step that
      ! was taken (rather then building MU based on step size in DMU)
      IC=2.
      !write(*,*) "IC:"
      !write(*,*) IC
      DO L=48,NPA-1 ! Fill the rest of the arrays not already filled
        PA(L+1)=CONE(IC) ! Start filling PA with values taken from the loss
        ! cone array
        IF(L.EQ.49) THEN
          PA(50)=16. ! hard code this angle as well
        ELSE
          ! This IFELSE block will slowly push IC to a value of SOMETHING, it
          ! always crashes when IC is 25 (while NR = 20). ALSO this if and
          ! else statement are behaving the same (i.e. IC will always increment
          ! by 1).
          if (IC.lt.nR) then
             !write(*,*) "we went to if"
             IC=IC+(nR-1.)/19.
          else
             !write(*,*) "we went to else"
             IC=IC+1
          endif
          !write(*,*) "IC:"
          !write(*,*) IC
        ENDIF
        MU(L+1)=COSD(PA(L+1)) ! define MU based on PA, filled from loss cone
        DMU(L)=(MU(L+1)-MU(L))       ! Grid size in cos pitch angle again back fill
        WMU(L)=2.*(DMU(L-1)-0.5*WMU(L-1))
        IF (L.GT.55) WMU(L)=0.5*(DMU(L)+DMU(L-1)) ! 55 is the value of Slen
      END DO
      DMU(NPA)=DMU(NPA-1)
      WMU(NPA)=DMU(NPA-1)
      ! When RadiusMax is not 6.5 the following DO loop will cause
      ! a crash durning the initilization step as it will eventually pass
      ! a value to acosd greater then 1. Currently the loop will not cause a
      ! crash if 5.5 < RadiusMax < 7.0. The crash does not happen after a
      ! consistent number of calls to acosd, sometimes only ~10 calls are made
      ! before passing an invaild value, others ~40 calls occur before a crash.
      DO L=1,NPA-1
        MUBOUN=MU(L)+0.5*WMU(L)
        if (MUBOUN.gt.1) then
          MUBOUN = 1
        else
          continue
        endif
        PAbn(L)=ACOSD(MUBOUN) ! PA at boundary of grid
      ENDDO
      PAbn(NPA)=0.
      write(*,*) "Print Statements form SWF_FILES make test1 of value of variables in src/ModRamInit.ARRAYS (f90)"
      write(*,*) "DL1"
      write(*,*) DL1
      write(*,*) "amla"
      write(*,*) amla
      write(*,*) "MDR"
      write(*,*) MDR
      write(*,*) "camlra"
      write(*,*) camlra
      write(*,*) "BE"
      write(*,*) BE
      write(*,*) "LZ, the wise"
      write(*,*) LZ
      write(*,*) "RLZ, the mystierous"
      write(*,*) RLZ
      write(*,*) "DPHI"
      write(*,*) DPHI
      write(*,*) "PHI"
      write(*,*) PHI
      write(*,*) "MLT"
      write(*,*) MLT
      write(*,*) "ELB"
      write(*,*) ELB
      write(*,*) "WE"
      write(*,*) WE
      write(*,*) "RW"
      write(*,*) RW
      write(*,*) "EKEV"
      write(*,*) EKEV
      write(*,*) "GREL"
      write(*,*) GREL
      write(*,*) "V"
      write(*,*) V
      write(*,*) "EBND"
      write(*,*) EBND
      write(*,*) "GRBND"
      write(*,*) GRBND
      write(*,*) "VBND"
      write(*,*) VBND
      write(*,*) "DE"
      write(*,*) DE
      write(*,*) "CONE, the dumb"
      write(*,*) CONE
      write(*,*) "PA"
      write(*,*) PA
      write(*,*) "MU"
      write(*,*) MU
      write(*,*) "WMU"
      write(*,*) WMU
      write(*,*) "DMU"
      write(*,*) DMU
      write(*,*) "PAbn"
      write(*,*) PAbn
    endif

    ! Determine the range of NPA such that PA is outside the loss cone:
    ! UPA is upper boundary for pitch angle for given Z
    ! When updated, UPA will need to be a function of radius and MLT -ME
     DO I=1,NR
       UPA(I) = NPA ! SZ, otherwise UPA = 0 for small enough loss cones
       DO L=NPA,1,-1
         IF(PA(L).LE.CONE(I)) UPA(I) = L     ! F(UPA)=0. - in loss cone
       END DO
     END DO

    ! calculate pitch angles for mlat
    DO I=1,NR
       DO IML=1,Slen
          DO IP=1,NPA
             spa=SQRT(SIND(PAbn(ip))**2*BE(i,iml)/BE(i,1))
             IF (spa.GT.1.0) spa=1.0
             ZRpabn(i,ip,iml)=ASIN(spa)
             IF (abs(spa-1.0).le.1e-9) THEN
                ZRpabn(i,ip,iml)=-1.0
             END IF
          ENDDO
       ENDDO
    ENDDO

    ! FFACTOR is ratio of F2 in conservative space to flux
    ! E* are factors to calculate temperature anisotropy
    DO I=1,NR
      DO K=1,NE
        DO L=2,NPA
          FFACTOR(S,I,K,L)=LZ(I)*LZ(I)*GREL(S,K)/SQRT(GREL(S,K)**2-1.)*MU(L)
          if (ffactor(s,i,k,l).le.0) print*,'s,i,k,l,ffactor=',s,i,k,l,ffactor(s,i,k,l)
        ENDDO
        FFACTOR(S,I,K,1)=FFACTOR(S,I,K,2)
      END DO
    END DO

    ! Energy factors used in RAM equations
    DO K=1,NE
      ERNH(S,K)=WE(K)*GREL(S,K)/SQRT((GREL(S,K)-1.)*(GREL(S,K)+1.)) ! [1/cm3]
      EPP(S,K)=ERNH(S,K)*EKEV(K)
      FACGR(S,K)=GREL(S,K)*SQRT((GREL(S,K)-1.)*(GREL(S,K)+1.))
    END DO

    ! to keep F constant at boundary
    CONF1=((LZ(NR)+DL1)/LZ(NR))**2
    CONF2=((LZ(NR)+2.*DL1)/LZ(NR))**2

    RFACTOR=3.4027E10*MDR*DPHI

    DEALLOCATE(CONE,RLAMBDA)
    RETURN
  END SUBROUTINE ARRAYS

!==================================================================================================
  SUBROUTINE init_input
    !!!! Module Variables
    use ModRamMain,      ONLY: DP, nIter
    use ModRamParams,    ONLY: IsRestart, IsStarttimeSet, NameBoundMag, &
                               DoUsePlasmasphere, HardRestart, InitializeOnFile, &
                               PlasmasphereModel
    use ModRamGrids,     ONLY: NL, NLT, nR, nT, nS
    use ModRamTiming,    ONLY: DtEfi, TimeRamNow, TimeRamElapsed
    use ModRamVariables, ONLY: Kp, F107, AE, TOLV, NECR, IP1, IR1, XNE, F2, species
    use ModScbParams,    ONLY: method, constTheta
    !!!! Module Subroutines/Functions
    use ModRamRun,       ONLY: ANISCH
    use ModRamIO,        ONLY: write_prefix
    use ModRamBoundary,  ONLY: get_boundary_flux
    use ModRamRestart,   ONLY: read_restart
    use ModRamIndices,   ONLY: get_indices
    use ModRamIO,        ONLY: read_initial, write2DFlux, writeLosses
    use ModRamInjection, ONLY: load_injection_file
    use ModRamFunctions, ONLY: ram_sum_pressure
    use ModRamScb,       ONLY: computehI, compute3DFlux
    use ModScbRun,       ONLY: scb_run, pressure
    use ModScbEuler,     ONLY: psiges, alfges
    use ModScbIO,        ONLY: computational_domain
    use ModScbCompute,   ONLY: computeBandJacob, compute_convergence
    use ModRamPlasmasphere, ONLY: plasmasphere
    !!!! Share Modules
    use ModIOUnit,      ONLY: UNITTMP_
    use ModTimeConvert, ONLY: TimeType

    implicit none

    integer :: i, j, j1, i1, iS, methodTemp

    character(len=4)   :: NameBoundMagTemp
    character(len=100) :: pmt
    character(len=100) :: HEADER
    character(len=*), parameter :: NameSub='init_input'


    !!!!!!!!!! Restart vs Initial Run
    if(IsRestart) then
       ! If Restart, read restart params and set timings appropriately.
       if (IsStarttimeSet) call CON_stop(NameSub//&
            ': Cannot use #STARTTIME command with #RESTART!')

       !!!!!! RESTART DATA !!!!!!!
       call read_restart

       call psiges
       call alfges

       call get_indices(TimeRamNow%Time, Kp, f107, AE)
       TOLV = FLOOR(TimeRamElapsed/DtEfi)*DtEfi

       ! Compute information not stored in restart files
       if (HardRestart) then
          call computational_domain
          call ram_sum_pressure
          call scb_run(0)
          call computehI(0)
          call compute3DFlux
       else
          call computeBandJacob
          call compute3DFlux
       endif

       call get_boundary_flux ! FGEOS
    else
       nIter = 1

       !!!!!! INITIALIZE DATA !!!!!
       ! Initial indices
       call get_indices(TimeRamNow%Time, Kp, f107, AE)
       TOLV = 0.0

       ! Initialize flux and pressure for default species
       ! Electron, Hydrogen, HeliumP1, OxygenP1 (sets to
       ! 0._dp for other species)
       call read_initial

       ! Initialize Magnetic Field (needed for injection files)
       ! ! If running with SWMF initialize on a dipole field
       if (NameBoundMag == 'SWMF') then
          NameBoundMagTemp = NameBoundMag
          methodTemp = method
          NameBoundMag = 'DIPL'
          method = 3
       endif
       call computational_domain

       ! Compute the SCB computational domain
       write(*,*) ''
       call write_prefix
       write(*,'(a)') 'Running SCB model to initialize B-field...'
       call ram_sum_pressure
       call scb_run(0)

       ! Couple SCB -> RAM
       call computehI(0)

       ! Initialize flux for species with injection files
       write(*,*) "About to initialize injections! punk."
       do i = 1, nS
          write(*,*) "i, Initfile = ", i, species(i)%Initialization
          if ((trim(species(i)%Initialization).ne.'InitializationFile') &
              .and.(trim(species(i)%Initialization).ne.'na')) then
             write(*,*) species(i)%Initialization
             call load_injection_file(i, trim(species(i)%Initialization), F2(i,:,:,:,:))
          endif
       enddo

       ! Reset dipole initialization for SWMF runs
       if (NameBoundMagTemp == 'SWMF') then
        method = methodTemp
        NameBoundMag = NameBoundMagTemp
       endif

       ! Plasmasphere Initialization
       if (DoUsePlasmasphere) then
          ! Read ig_rz.dat and apf107.dat
          call read_ig_rz
          call readapf107

          pmt = PlasmasphereModel
          PlasmasphereModel = "Carpenter"
          call plasmasphere(0._dp)
          PlasmasphereModel = pmt
       endif

       ! Call initial outputs
       call compute3DFlux
       call write2DFlux
       call writeLosses
    end if
  !!!!!!!!

   return

  end subroutine init_input

END MODULE ModRamInit
