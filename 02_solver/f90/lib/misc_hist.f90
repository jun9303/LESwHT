!=======================================================================
       subroutine ftrfiles(io)
!=======================================================================
!$       USE OMP_LIB
         use mod_common
         implicit none
         integer(8) :: io

         if (io .ne. 0) then
           if (iread .ne. 1) then
             open (2000, file='../output/ftr/fhist.dat')
             open (2001, file='../output/ftr/fcdcl.dat')
             open (2002, file='../output/ftr/fnusselt.dat')
             if (ntrace .ne. 0) then
               open (2005, file='../output/ftr/futrace.dat')
               open (2006, file='../output/ftr/fvtrace.dat')
               open (2007, file='../output/ftr/fwtrace.dat')
               open (2008, file='../output/ftr/fptrace.dat')
             end if
             if (ich .ne. 0) then
               open (2010, file='../output/ftr/fcmfr.dat')
             end if
             open (2011, file='../output/ftr/ftime.dat')
           else
             open (2000, file='../output/ftr/fhist.dat', &
                   position='append')
             open (2001, file='../output/ftr/fcdcl.dat', &
                   position='append')
             open (2002, file='../output/ftr/fnusselt.dat', &
                   position='append')
             if (ntrace .ne. 0) then
               open (2005, file='../output/ftr/futrace.dat', &
                     position='append')
               open (2006, file='../output/ftr/fvtrace.dat', &
                     position='append')
               open (2007, file='../output/ftr/fwtrace.dat', &
                     position='append')
               open (2008, file='../output/ftr/fptrace.dat', &
                     position='append')
             end if
             if (ich .ne. 0) then
               open (2010, file='../output/ftr/fcmfr.dat', &
                     position='append')
             end if
             open (2011, file='../output/ftr/ftime.dat', &
                   position='append')
           end if
         else
           close (2000)
           close (2001)
           close (2002)
           if (ntrace .ne. 0) then
             close (2005)
             close (2006)
             close (2007)
             close (2008)
           end if
           if (ich .ne. 0) then
             close (2010)
           end if
           close (2011)

         end if

         return
       end subroutine ftrfiles
!=======================================================================
!=======================================================================
       subroutine writehistory
!=======================================================================
!
!     TRACE THE MAXIMUM CFL#, DIVERGENCEU, AND QMASS(IBM).
!     THE TRACING INTERVAL IS SET IN LICA.IN (NPIN)
!
!     *IN CASE OF CHANNEL FLOW,
!     - THE FLOW RATE AND THE MEAN PRESSURE GRADIENT IS TRACED
!     - QFLUX SHOULD BE CONSTANT, WHILE PMI, PMIL, PMIU FLUCTUATE
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray
         implicit none
         integer(8) :: i, j, k
         real(8) :: ubulk_inst, cmfravg(3), funcbody, flowvol, qvol_u

         qflux = 0.
         flowvol = 0.0d0
         qvol_u = 0.0d0

!$OMP PARALLEL DO &
!$OMP REDUCTION(+:QVOL_U,FLOWVOL)
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m
               if (funcbody(x(i), ymp(j), zmp(k), time) .ge. 1.e-10) then
                 qvol_u = qvol_u + u(i, j, k) * c2cx(i) * f2fy(j) * f2fz(k)
                 flowvol = flowvol + c2cx(i) * f2fy(j) * f2fz(k)
               end if
             end do
           end do
         end do

         qflux = qvol_u / xl
         if (flowvol .gt. 0.0d0) then
           ubulk_inst = qvol_u / flowvol
         else
           ubulk_inst = 0.0d0
         end if
         cmfravg = 0.0d0

         write (2000, 130) time, dt, cflmax, dvmax, qmmax
130      format(f13.5, 4es15.7)

         if (ich .ne. 0) then
           if (time .ge. avg_tst) then
             cmfravg_int(1) = cmfravg_int(1) + qflux * dt
             cmfravg_int(2) = cmfravg_int(2) + ubulk_inst * dt
             cmfravg_int(3) = cmfravg_int(3) + pmiavg * dt
             cmfravg_dur = cmfravg_dur + dt
             if (cmfravg_dur .gt. 0.0d0) then
               cmfravg(1) = cmfravg_int(1) / cmfravg_dur
               cmfravg(2) = cmfravg_int(2) / cmfravg_dur
               cmfravg(3) = cmfravg_int(3) / cmfravg_dur
             end if
           end if

           write (2010, 140) time, qflux, ubulk_inst, pmiavg, cmfravg(1), cmfravg(2), cmfravg(3)
140        format(f13.5, 6es15.7)
         end if

         return
       end subroutine writehistory
!=======================================================================
!=======================================================================
       subroutine tracer
!=======================================================================
!
!     THIS SUBROUTINE IS CALLED WHEN MOD(M,NTPRINT)=0.
!     VELOCITY COMPONENTS AND THE PRESSURE AT THE CELL CENTER OF THE
!     TRACING CELLS ARE SAVED.
!
!     NTPRINT : DETERMINING THE TIME-TRACING INTERVAL
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: u, v, w, p
         implicit none
         integer(8) :: n

         write (2005, 101) time, (0.5 * (u(trpts(n, 1), trpts(n, 2), trpts(n, 3)) &
                                         + u(trpts(n, 1) + 1, trpts(n, 2), trpts(n, 3))), n=1, ntrace)
         write (2006, 101) time, (0.5 * (v(trpts(n, 1), trpts(n, 2), trpts(n, 3)) &
                                         + v(trpts(n, 1), trpts(n, 2) + 1, trpts(n, 3))), n=1, ntrace)
         write (2007, 101) time, (0.5 * (w(trpts(n, 1), trpts(n, 2), trpts(n, 3)) &
                                         + w(trpts(n, 1), trpts(n, 2), trpts(n, 3) + 1)), n=1, ntrace)
         write (2008, 101) time, (p(trpts(n, 1), trpts(n, 2), trpts(n, 3)), n=1, ntrace)
101      format(f15.7, 10000es15.7)

         return
       end subroutine tracer
!=======================================================================
!=======================================================================
       subroutine writefield
!=======================================================================
!
!     MAKE AN OUTPUT FILE OF INSTANTANEOUS FIELD, WHEN MOD(NTIME,NPRINT)=0,
!     NPRINT : INSTANTANEOUS FIELD FILE PRINTING INTERVAL
!     TFN1='FLD': PREFIX FOR INSTANTANEOUS FLOW FIELD
!     TNAME     : INSTANTANEOUS FIELD FILE NAME EX) FLD0006100
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: u, v, w, t, p
         implicit none
         integer(8) :: i, j, k
         integer(8) :: idum
         real(8) :: dum
         character*26 :: tname
         character*19 :: tfn1

         idum = 0
         dum = 0.

         tfn1 = '../output/field/fld'
         write (tname, '(a,i7.7)') tfn1, ihist

         open (nv, file=tname, form='unformatted')
         write (nv) n1, n2, n3, re, pr, gr
         write (nv) ihist, m, time, dt
         write (nv) xprdic, yprdic, zprdic
         write (nv) bc_xbtm, bc_xtop, bc_ybtm, bc_ytop, bc_zbtm, bc_ztop
         write (nv) ich, iconjg
         write (nv) (((u(i, j, k), i=1, n1), j=0, n2), k=0, n3)
         write (nv) (((v(i, j, k), i=0, n1), j=1, n2), k=0, n3)
         write (nv) (((w(i, j, k), i=0, n1), j=0, n2), k=1, n3)
         write (nv) (((p(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
         if (ihtrans .eq. 1) write (nv) (((t(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
         close (nv)

         nv = nv + 1

         return
       end subroutine writefield
!=======================================================================
!=======================================================================
       subroutine writeftrtime
!=======================================================================
!
!     REAL TIME MEASUREMENTS FOR INDIVIDUAL PROCESSES
!     PRINTED AT EVERY COMPUTATIONAL STEP
!
!-----------------------------------------------------------------------
         use mod_common
         implicit none
         real(8) :: ftrtime1, ftrtime2, ftrtime3, ftrtime4

         ftrtime1 = (sgstime_e(3) - sgstime_b(3) + sgstime_e(2) - sgstime_b(2) + sgstime_e(1) - sgstime_b(1))
         ftrtime2 = (rhsnlhstime_e(3) - rhsnlhstime_b(3) + rhsnlhstime_e(2) - rhsnlhstime_b(2) + rhsnlhstime_e(1) - rhsnlhstime_b(1))
         ftrtime3 = (poisstime_e(3) - poisstime_b(3) + poisstime_e(2) - poisstime_b(2) + poisstime_e(1) - poisstime_b(1))
         ftrtime4 = time_end - time_begin

         write (*, 201) ftrtime1
         write (*, 202) ftrtime2
         write (*, 203) ftrtime3
         write (*, 204) ftrtime4

201      format('TIME FOR SGS    : ', f12.3, ' SECONDS')
202      format('TIME FOR RHSNLHS: ', f12.3, ' SECONDS')
203      format('TIME FOR POISSON: ', f12.3, ' SECONDS')
204      format('TIME OF OPERTION: ', f12.3, ' SECONDS')
         write (2011, 206) time, ftrtime1, ftrtime2, ftrtime3, ftrtime4
206      format(f13.5, 4es15.7)

         return
       end
!=======================================================================
!=======================================================================
       subroutine field_avg
!=======================================================================
!     SUBROUTINE FOR FIELD AVERAGING
!     AVERAGED VARIABLES ARE DEFINED AT CELL CENTER
!     VARIABLES ARE LINEARLY INTERPOLATED
!     FROM STAGGERED GRID STRUCTURE TO DETERMINE CELL CENTER VALUES
!
!     NO SPATIAL AVERAGING.
!
!     PRINT AN AVERAGED FLOW-FIELD FILE, WHEN MOD(NTIME,NPRIAVG)=0,
!     SATISFYING THE AVERAGED FLOW-FIELD FILE PRINTING INTERVAL (NPRIAVG SET IN LICA.IN).
!     TFN1='FAV': PREFIX FOR AVERAGE FLOW FIELD
!     TNAME     : AVERAGE FIELD FILE NAME  EX) FAV0100000-0110000
!
!-----------------------------------------------------------------------
         use mod_common
         use mod_flowarray, only: u, v, w, p, t, uavg, vavg, wavg, uiujavg, pavg, &
                                  p2avg, tavg, t2avg, voravg, vor2avg, ssavg
         implicit none
         integer(8) :: i, j, k, l
         real(8) :: ucc, vcc, wcc
         real(8) :: vg11, vg12, vg13, vg21, vg22, vg23, vg31, vg32, vg33
         real(8) :: up, um, vp, vm, wp, wm, sr(6), srsr
         real(8) :: wx, wy, wz
         real(8) :: timeend
         integer(8) :: ihistend
         character*18 :: tname
         character*3 :: tfn1
         character*7 :: tfn2, tfn3
         character*1 :: tfnh

         tfnh = '-'

         if (time .le. avg_tst) return

         if (.not. avg_started) then
           avg_started = .true.
           npriavg_count = 0
           ihistavg_start = ihist
           timeinit = avg_tst
           ihistinit = ihist
           uavg = 0.
           vavg = 0.
           wavg = 0.
           uiujavg = 0.
           pavg = 0.
           p2avg = 0.
           if (ihtrans .eq. 1) then
             tavg = 0.
             t2avg = 0.
           end if
           voravg = 0.
           vor2avg = 0.
           ssavg = 0.
         end if

!$OMP PARALLEL DO &
!$OMP PRIVATE(UCC,VCC,WCC)
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m
               ucc = 0.5 * (u(i, j, k) + u(i + 1, j, k))
               vcc = 0.5 * (v(i, j, k) + v(i, j + 1, k))
               wcc = 0.5 * (w(i, j, k) + w(i, j, k + 1))
               uavg(i, j, k) = uavg(i, j, k) + ucc * dt
               vavg(i, j, k) = vavg(i, j, k) + vcc * dt
               wavg(i, j, k) = wavg(i, j, k) + wcc * dt
               uiujavg(i, j, k, 1) = uiujavg(i, j, k, 1) + ucc**2.*dt
               uiujavg(i, j, k, 2) = uiujavg(i, j, k, 2) + ucc * vcc * dt
               uiujavg(i, j, k, 3) = uiujavg(i, j, k, 3) + ucc * wcc * dt
               uiujavg(i, j, k, 4) = uiujavg(i, j, k, 4) + vcc**2.*dt
               uiujavg(i, j, k, 5) = uiujavg(i, j, k, 5) + vcc * wcc * dt
               uiujavg(i, j, k, 6) = uiujavg(i, j, k, 6) + wcc**2.*dt
               if (ihtrans .eq. 1) then
                 tavg(i, j, k) = tavg(i, j, k) + t(i, j, k) * dt
                 t2avg(i, j, k) = t2avg(i, j, k) + t(i, j, k)**2.*dt
               end if
               pavg(i, j, k) = pavg(i, j, k) + p(i, j, k) * dt
               p2avg(i, j, k) = p2avg(i, j, k) + p(i, j, k)**2.*dt
             end do
           end do
         end do

!$OMP PARALLEL DO &
!$OMP PRIVATE(VG11,VG12,VG13,VG21,VG22,VG23,VG31,VG32,VG33) &
!$OMP PRIVATE(UP,UM,VP,VM,WP,WM,SR,SRSR) &
!$OMP PRIVATE(WX,WY,WZ)
         do k = 1, n3m
           do j = 1, n2m
             do i = 1, n1m

               vg11 = f2fxi(i) * (u(i + 1, j, k) - u(i, j, k))
               up = c2cyi(j + 1) * 0.25 &
                    * (f2fy(j + 1) * (u(i, j, k) + u(i + 1, j, k)) &
                       + f2fy(j) * (u(i, j + 1, k) + u(i + 1, j + 1, k))) * (1.-fixju(j)) &
                    + 0.5 * (u(i, n2, k) + u(i + 1, n2, k)) * fixju(j)
               um = c2cyi(j) * 0.25 &
                    * (f2fy(j) * (u(i, j - 1, k) + u(i + 1, j - 1, k)) &
                       + f2fy(j - 1) * (u(i, j, k) + u(i + 1, j, k))) * (1.-fixjl(j)) &
                    + 0.5 * (u(i, 0, k) + u(i + 1, 0, k)) * fixjl(j)
               vg12 = f2fyi(j) * (up - um)

               up = c2czi(k + 1) * 0.25 &
                    * (f2fz(k + 1) * (u(i, j, k) + u(i + 1, j, k)) &
                       + f2fz(k) * (u(i, j, k + 1) + u(i + 1, j, k + 1))) * (1.-fixku(k)) &
                    + 0.5 * (u(i, j, n3) + u(i + 1, j, n3)) * fixku(k)
               um = c2czi(k) * 0.25 &
                    * (f2fz(k) * (u(i, j, k - 1) + u(i + 1, j, k - 1)) &
                       + f2fz(k - 1) * (u(i, j, k) + u(i + 1, j, k))) * (1.-fixkl(k)) &
                    + 0.5 * (u(i, j, 0) + u(i + 1, j, 0)) * fixkl(k)
               vg13 = f2fzi(k) * (up - um)

               vp = c2cxi(i + 1) * 0.25 &
                    * (f2fx(i + 1) * (v(i, j, k) + v(i, j + 1, k)) &
                       + f2fx(i) * (v(i + 1, j, k) + v(i + 1, j + 1, k))) * (1.-fixiu(i)) &
                    + 0.5 * (v(n1, j, k) + v(n1, j + 1, k)) * fixiu(i)
               vm = c2cxi(i) * 0.25 &
                    * (f2fx(i) * (v(i - 1, j, k) + v(i - 1, j + 1, k)) &
                       + f2fx(i - 1) * (v(i, j, k) + v(i, j + 1, k))) * (1.-fixil(i)) &
                    + 0.5 * (v(0, j, k) + v(0, j + 1, k)) * fixil(i)
               vg21 = f2fxi(i) * (vp - vm)

               vg22 = f2fyi(j) * (v(i, j + 1, k) - v(i, j, k))

               vp = c2czi(k + 1) * 0.25 &
                    * (f2fz(k + 1) * (v(i, j, k) + v(i, j + 1, k)) &
                       + f2fz(k) * (v(i, j, k + 1) + v(i, j + 1, k + 1))) * (1.-fixku(k)) &
                    + 0.5 * (v(i, j, n3) + v(i, j + 1, n3)) * fixku(k)
               vm = c2czi(k) * 0.25 &
                    * (f2fz(k) * (v(i, j, k - 1) + v(i, j + 1, k - 1)) &
                       + f2fz(k - 1) * (v(i, j, k) + v(i, j + 1, k))) * (1.-fixkl(k)) &
                    + 0.5 * (v(i, j, 0) + v(i, j + 1, 0)) * fixkl(k)
               vg23 = f2fzi(k) * (vp - vm)

               wp = c2cxi(i + 1) * 0.25 &
                    * (f2fx(i + 1) * (w(i, j, k) + w(i, j, k + 1)) &
                       + f2fx(i) * (w(i + 1, j, k) + w(i + 1, j, k + 1))) * (1.-fixiu(i)) &
                    + 0.5 * (w(n1, j, k) + w(n1, j, k + 1)) * fixiu(i)
               wm = c2cxi(i) * 0.25 &
                    * (f2fx(i) * (w(i - 1, j, k) + w(i - 1, j, k + 1)) &
                       + f2fx(i - 1) * (w(i, j, k) + w(i, j, k + 1))) * (1.-fixil(i)) &
                    + 0.5 * (w(0, j, k) + w(0, j, k + 1)) * fixil(i)
               vg31 = f2fxi(i) * (wp - wm)

               wp = c2cyi(j + 1) * 0.25 &
                    * (f2fy(j + 1) * (w(i, j, k) + w(i, j, k + 1)) &
                       + f2fy(j) * (w(i, j + 1, k) + w(i, j + 1, k + 1))) * (1.-fixju(j)) &
                    + 0.5 * (w(i, n2, k) + w(i, n2, k + 1)) * fixju(j)
               wm = c2cyi(j) * 0.25 &
                    * (f2fy(j) * (w(i, j - 1, k) + w(i, j - 1, k + 1)) &
                       + f2fy(j - 1) * (w(i, j, k) + w(i, j, k + 1))) * (1.-fixjl(j)) &
                    + 0.5 * (w(i, 0, k) + w(i, 0, k + 1)) * fixjl(j)
               vg32 = f2fyi(j) * (wp - wm)

               vg33 = f2fzi(k) * (w(i, j, k + 1) - w(i, j, k))

               sr(1) = vg11
               sr(2) = 0.5 * (vg12 + vg21)
               sr(3) = 0.5 * (vg13 + vg31)
               sr(4) = vg22
               sr(5) = 0.5 * (vg23 + vg32)
               sr(6) = vg33

               wx = vg32 - vg23
               wy = vg13 - vg31
               wz = vg21 - vg12

               voravg(i, j, k, 1) = voravg(i, j, k, 1) + wx * dt
               voravg(i, j, k, 2) = voravg(i, j, k, 2) + wy * dt
               voravg(i, j, k, 3) = voravg(i, j, k, 3) + wz * dt

               vor2avg(i, j, k, 1) = vor2avg(i, j, k, 1) + wx * wx * dt
               vor2avg(i, j, k, 2) = vor2avg(i, j, k, 2) + wx * wy * dt
               vor2avg(i, j, k, 3) = vor2avg(i, j, k, 3) + wx * wz * dt
               vor2avg(i, j, k, 4) = vor2avg(i, j, k, 4) + wy * wy * dt
               vor2avg(i, j, k, 5) = vor2avg(i, j, k, 5) + wy * wz * dt
               vor2avg(i, j, k, 6) = vor2avg(i, j, k, 6) + wz * wz * dt

               srsr = 2.*sr(2)**2.+sr(1)**2.+2.*sr(3)**2.+sr(4)**2.+2.*sr(5)**2.+sr(6)**2.
               ssavg(i, j, k) = ssavg(i, j, k) + srsr * dt

             end do
           end do
         end do

         npriavg_count = npriavg_count + 1

         if (npriavg_count .ge. npriavg) then

           timeend = time
           ihistend = ihist

           tfn1 = 'fav'
           write (tfn2, '(i7.7)') ihistinit
           write (tfn3, '(i7.7)') ihist
           tname = tfn1//tfn2//tfnh//tfn3

           open (nav, file='../output/field_avg/'//tname, form='unformatted')
           write (nav) n1m, n2m, n3m, re
          write (nav) timeinit, timeend, dble(ihistinit), dble(ihistend)
           write (nav) (((uavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
           write (nav) (((vavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
           write (nav) (((wavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
           write (nav) ((((uiujavg(i, j, k, l), i=1, n1m), j=1, n2m), k=1, n3m), l=1, 6)
           write (nav) (((pavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
           write (nav) (((p2avg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
           write (nav) ((((voravg(i, j, k, l), i=1, n1m), j=1, n2m), k=1, n3m), l=1, 3)
           write (nav) ((((vor2avg(i, j, k, l), i=1, n1m), j=1, n2m), k=1, n3m), l=1, 6)
           write (nav) (((ssavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
           if (ihtrans .eq. 1) then
             write (nav) (((tavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
             write (nav) (((t2avg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
           end if
           close (nav)

           open (2999, file='../output/field_avg/fav_manifest.dat', position='append')
           write (2999, '(a)') trim(tname)
           close (2999)

           nav = nav + 1
           uavg = 0.
           vavg = 0.
           wavg = 0.
           uiujavg = 0.
           pavg = 0.
           p2avg = 0.
           if (ihtrans .eq. 1) then
             tavg = 0.
             t2avg = 0.
           end if
           voravg = 0.
           vor2avg = 0.
           ssavg = 0.
           npriavg_count = 0
           timeinit = time
           ihistinit = ihist

         end if

         return
       end
!=======================================================================
!=======================================================================
      subroutine field_avg_finalize
!=======================================================================
        use mod_common
        use mod_flowarray, only: uavg, vavg, wavg, uiujavg, pavg, p2avg, tavg, t2avg, voravg, vor2avg, ssavg
        implicit none
        integer(8) :: i, j, k, l, iu, ios, file_count
        integer :: cmdstat
        integer(8) :: n1mr, n2mr, n3mr, ihist0, ihist1
        integer(8) :: ihistend
        real(8) :: rer, time0, time1, timeend_all
        logical :: have_data
        character*24 :: tname
        character*24 :: fname
        character*7 :: tfn1
        character*7 :: tfn2, tfn3
        character*1 :: tfnh
        real(8), allocatable :: tmp3d(:, :, :), tmp4d3(:, :, :, :), tmp4d6(:, :, :, :)

        if ((iavg .ne. 1) .or. (.not. avg_started)) return

        tfnh = '-'

        if (npriavg_count .gt. 0) then
          ihistend = ihist

          write (tfn2, '(i7.7)') ihistinit
          write (tfn3, '(i7.7)') ihist
          tname = 'fav'//tfn2//tfnh//tfn3

          open (nav, file='../output/field_avg/'//tname, form='unformatted')
          write (nav) n1m, n2m, n3m, re
          write (nav) timeinit, time, dble(ihistinit), dble(ihistend)
          write (nav) (((uavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
          write (nav) (((vavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
          write (nav) (((wavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
          write (nav) ((((uiujavg(i, j, k, l), i=1, n1m), j=1, n2m), k=1, n3m), l=1, 6)
          write (nav) (((pavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
          write (nav) (((p2avg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
          write (nav) ((((voravg(i, j, k, l), i=1, n1m), j=1, n2m), k=1, n3m), l=1, 3)
          write (nav) ((((vor2avg(i, j, k, l), i=1, n1m), j=1, n2m), k=1, n3m), l=1, 6)
          write (nav) (((ssavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
          if (ihtrans .eq. 1) then
            write (nav) (((tavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
            write (nav) (((t2avg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
          end if
          close (nav)

          open (2999, file='../output/field_avg/fav_manifest.dat', position='append')
          write (2999, '(a)') trim(tname)
          close (2999)
          nav = nav + 1
        end if

        uavg = 0.
        vavg = 0.
        wavg = 0.
        uiujavg = 0.
        pavg = 0.
        p2avg = 0.
        if (ihtrans .eq. 1) then
          tavg = 0.
          t2avg = 0.
        end if
        voravg = 0.
        vor2avg = 0.
        ssavg = 0.

        allocate (tmp3d(n1m, n2m, n3m))
        allocate (tmp4d3(n1m, n2m, n3m, 3))
        allocate (tmp4d6(n1m, n2m, n3m, 6))

        have_data = .false.
        file_count = 0
        timeend_all = avg_tst
        open (2998, file='../output/field_avg/fav_manifest.dat', status='old', iostat=ios)
        if (ios .ne. 0) then
          deallocate (tmp3d, tmp4d3, tmp4d6)
          return
        end if

        do
          read (2998, '(a)', iostat=ios) fname
          if (ios .ne. 0) exit
          if (len_trim(fname) .le. 0) cycle

          iu = 2997
          open (iu, file='../output/field_avg/'//trim(fname), form='unformatted', status='old', iostat=ios)
          if (ios .ne. 0) cycle

          read (iu, iostat=ios) n1mr, n2mr, n3mr, rer
          if (ios .ne. 0) then
            close (iu)
            cycle
          end if
          read (iu) time0, time1, ihist0, ihist1

          read (iu) tmp3d
          uavg = uavg + tmp3d
          read (iu) tmp3d
          vavg = vavg + tmp3d
          read (iu) tmp3d
          wavg = wavg + tmp3d
          read (iu) tmp4d6
          uiujavg = uiujavg + tmp4d6
          read (iu) tmp3d
          pavg = pavg + tmp3d
          read (iu) tmp3d
          p2avg = p2avg + tmp3d
          read (iu) tmp4d3
          voravg = voravg + tmp4d3
          read (iu) tmp4d6
          vor2avg = vor2avg + tmp4d6
          read (iu) tmp3d
          ssavg = ssavg + tmp3d
          if (ihtrans .eq. 1) then
            read (iu) tmp3d
            tavg = tavg + tmp3d
            read (iu) tmp3d
            t2avg = t2avg + tmp3d
          end if
          close (iu)

          file_count = file_count + 1
          if (.not. have_data) then
            timeinit = time0
            ihistinit = ihist0
            have_data = .true.
          end if
          timeend_all = time1
          ihistend = ihist1
        end do
        close (2998)

        deallocate (tmp3d, tmp4d3, tmp4d6)

        if ((.not. have_data) .or. (file_count .le. 0)) return

        write (tfn2, '(i7.7)') ihistinit
        write (tfn3, '(i7.7)') ihistend
        tfn1 = 'fav'
        tname = tfn1//tfn2//tfnh//tfn3

        call execute_command_line('mkdir -p ../output/field_avg/combined', exitstat=cmdstat)
        if (cmdstat .ne. 0) then
          write (*, *) 'WARNING: failed to create ../output/field_avg/combined; writing combined file in ../output/field_avg/'
          open (nav, file='../output/field_avg/'//tname, form='unformatted')
        else
          open (nav, file='../output/field_avg/combined/'//tname, form='unformatted')
        end if
        write (nav) n1m, n2m, n3m, re
        write (nav) timeinit, timeend_all, dble(ihistinit), dble(ihistend)
        write (nav) (((uavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
        write (nav) (((vavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
        write (nav) (((wavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
        write (nav) ((((uiujavg(i, j, k, l), i=1, n1m), j=1, n2m), k=1, n3m), l=1, 6)
        write (nav) (((pavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
        write (nav) (((p2avg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
        write (nav) ((((voravg(i, j, k, l), i=1, n1m), j=1, n2m), k=1, n3m), l=1, 3)
        write (nav) ((((vor2avg(i, j, k, l), i=1, n1m), j=1, n2m), k=1, n3m), l=1, 6)
        write (nav) (((ssavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
        if (ihtrans .eq. 1) then
          write (nav) (((tavg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
          write (nav) (((t2avg(i, j, k), i=1, n1m), j=1, n2m), k=1, n3m)
        end if
        close (nav)
        nav = nav + 1

        return
      end subroutine field_avg_finalize
!=======================================================================
