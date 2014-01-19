!                                                                                      
!  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”        
!  or “3-clause license”)                                                              
!  Please read attached file License.txt                                               
!                                        
      subroutine timer(ttime)
      use const
      real(dp) :: ttime
!
      real(sp) :: temp
!
!     This routine computes cpu time in double precision; it makes use of 
!     the intrinsic f90 cpu_time therefore a conversion type is
!     needed.
!
!           J.L Morales  Departamento de Matematicas, 
!                        Instituto Tecnologico Autonomo de Mexico
!                        Mexico D.F.
!
!           J.L Nocedal  Department of Electrical Engineering and
!                        Computer Science.
!                        Northwestern University. Evanston, IL. USA
!                         
!                        January 21, 2011
!
      temp = sngl(ttime)
      call cpu_time(temp)
      ttime = dble(temp) 

      return

      end
      
