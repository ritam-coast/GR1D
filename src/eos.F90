!-*-f90-*-
subroutine eos_full(i,xrhoi,    &
     xtemp,                   &
     xye,                     &
     xenri,                   &
     xprs,                    &
     xprs_th,                 & 
     xent,                    &
     xcs2,                    &
     xdedt,                   &
     xdpderho,                &
     xdpdrhoe,                & 
     xxa,xxh,xxn,xxp,         & 
     xabar,xzbar,             &
     xmu_e,xmu_n,xmu_p,xmunu,& 
     keytemp,keyerr,eoskey,rfeps)

  use GR1D_module, only: clite,rho_gf,press_gf,eps_gf,GR,n1,atmo
  use atmos
  implicit none

  integer,intent(in)    :: i
  real*8, intent(inout) :: xrhoi  ! inout because of atmosphere
  real*8, intent(in)    :: xye
  real*8, intent(inout) :: xtemp,xenri,xent
  real*8, intent(out)   :: xprs,xprs_th,xcs2,xdedt
  real*8, intent(out)   :: xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp
  real*8, intent(out)   :: xabar,xzbar,xmu_e,xmu_n,xmu_p,xmunu
  real*8, intent(in)    :: rfeps
  integer, intent(in)   :: keytemp
  integer, intent(in)   :: eoskey
  integer, intent(out)  :: keyerr

  real*8 :: xrho, xenr
  

  ! convert to CGS:
  xrho = xrhoi/rho_gf
  xenr = xenri/eps_gf

!  if(keyerr.eq.0) then
     if(xrhoi.le.atmo_rho_thr) then
        call atmos_eos(i,xrhoi,xprs,xenri,xcs2)
        return
     endif
!  endif

  if(eoskey.eq.1) then
     ! hybrid EOS
     call hybrid_eos(xrho,xenr,xprs,xprs_th,&
          xdpdrhoe,xdpderho,xcs2,keytemp)
     xprs = xprs*press_gf
     xprs_th = xprs_th*press_gf
     if (GR) then
        xcs2 = xcs2/(clite**2*(1.0d0+xprs/xrhoi+xenri))
     else
        xcs2 = xcs2/clite**2
     endif
     xdpderho = xdpderho*press_gf/eps_gf
     xdpdrhoe = xdpdrhoe*press_gf/rho_gf
     xent = 0.0d0
     xdedt = 0.0d0
     xxa = 0.0d0
     xxh = 0.0d0
     xxp = 0.0d0
     xxn = 0.0d0
     xabar = 0.0d0
     xzbar = 0.0d0
     xmu_e = 0.0d0
     xmu_n = 0.0d0
     xmu_p = 0.0d0
     xmunu = 0.0d0
     if(keytemp.eq.1) then
        xenri = xenr*eps_gf
     endif
     keyerr = 0
  else if(eoskey.eq.4) then
     ! Gamma-Law EOS
     call ideal_eos(xrho,xenr,xprs,xdpdrhoe,xdpderho,xcs2,keytemp)
     xprs = xprs*press_gf
     xprs_th = xprs
     if (GR) then
        xcs2 = xcs2/(clite**2*(1.0d0+xprs/xrhoi+xenri))
     else
        xcs2 = xcs2/clite**2
     endif
     xdpderho = xdpderho*press_gf/eps_gf
     xdpdrhoe = xdpdrhoe*press_gf/rho_gf
     xent = 0.0d0
     xdedt = 0.0d0
     xxa = 0.0d0
     xxh = 0.0d0
     xxp = 0.0d0
     xxn = 0.0d0
     xabar = 0.0d0
     xzbar = 0.0d0
     xmu_e = 0.0d0
     xmu_n = 0.0d0
     xmu_p = 0.0d0
     xmunu = 0.0d0
     if(keytemp.eq.1) then
        xenri = xenr*eps_gf
     endif
     keyerr = 0
  else if(eoskey.eq.2) then
     ! Poly EOS
     call poly_eos(xrho,xenr,xprs,xdpdrhoe,xdpderho,xcs2,keytemp)
     xprs = xprs*press_gf
     xprs_th = xprs
     if (GR) then
        xcs2 = xcs2/(clite**2*(1.0d0+xprs/xrhoi+xenri))
     else
        xcs2 = xcs2/clite**2
     endif
     xdpderho = xdpderho*press_gf/eps_gf
     xdpdrhoe = xdpdrhoe*press_gf/rho_gf
     xent = 0.0d0
     xdedt = 0.0d0
     xxa = 0.0d0
     xxh = 0.0d0
     xxp = 0.0d0
     xxn = 0.0d0
     xabar = 0.0d0
     xzbar = 0.0d0
     xmu_e = 0.0d0
     xmu_n = 0.0d0
     xmu_p = 0.0d0
     xmunu = 0.0d0
     if(keytemp.eq.1) then
        xenri = xenr*eps_gf
     endif
     keyerr = 0
  else if(eoskey.eq.3) then
#if HAVE_NUC_EOS
     call nuc_eos_short(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
          xdpderho,xdpdrhoe,xmunu,keytemp,keyerr,rfeps)
     xprs = xprs*press_gf
     xprs_th = 0.0d0
     if (GR) then
        xcs2 = xcs2/(clite**2*(1.0d0+xprs/xrhoi+xenri))
     else
        xcs2 = xcs2/clite**2
     endif
     xdpderho = xdpderho*press_gf/eps_gf
     xdpdrhoe = xdpdrhoe*press_gf/rho_gf
     xxa = 0.0d0
     xxh = 0.0d0
     xxp = 0.0d0
     xxn = 0.0d0
     xabar = 0.0d0
     xzbar = 0.0d0
     xmu_e = 0.0d0
     xmu_n = 0.0d0
     xmu_p = 0.0d0
     xmunu = xmunu
     if(keytemp.eq.1.or.keytemp.eq.2) then
        xenri = xenr*eps_gf
     endif
     if(keyerr.eq.667) then
        if(i.ne.n1) then
           if(atmo(i+1).ne.0) then
              atmo(i) = 1
              call atmos_eos(i,xrhoi,xprs,xenri,xcs2)
              keyerr = 0
              stop "eh"
           endif
        endif
     endif
#else
     stop "eoskey=3 impossible, since NUC_EOS not present."
#endif
  else
     write(6,*) "eoskey ",eoskey," not implemented with eos_full"
     stop "This is bad! Fix me please!"
  endif


end subroutine eos_full

subroutine eos(i,ri,tio,y,eio,xx,keytemp,keyerr,eosflag,eoskey,rfeps)

  use GR1D_module, only: clite,rho_gf,press_gf,eps_gf,GR,n1,atmo
  use atmos
  implicit none

  integer, intent(in) :: i
  real*8, intent(in) :: y,rfeps
  real*8, intent(inout) :: ri,eio, tio
  real*8, intent(out)    :: xx

  real*8 r,e,p_th,tp
  real*8 px,sx,cs2x,gammax,dedt,ex
  integer eosflag,eoskey,keytemp,keyerr,keyerrt
  !internal
  real*8 prs,soundsqr
  real*8 dpdrho,dpde

#if HAVE_NUC_EOS
  real*8 xmunu,xent,xdedt
#endif

  xx = 0.0d0

  ! convert to CGS:
  r = ri/rho_gf
  e = eio/eps_gf

! eosflags:                                                                     
! 1 --> pressure                                                                
! 2 --> dpdrhoe                                                                 
! 3 --> dpderho                                                                 
! 4 --> eps                                                                     
! 5 --> temp                                                                    
! 6 --> cs2                                                                     
! 7 --> gamma                                                                   
! 8 --> entropy                                                                 
! 9 --> munu             

! eoskey:
! 1 --> Hybrid EOS
! 2 --> Poly EOS
! 3 --> finite-T EOS
! 4 --> Gamma-Law EOS

     if(ri.le.atmo_rho_thr) then
        call atmos_eos(i,ri,prs,eio,soundsqr)
        
        select case (eosflag)
        case(1)
           !pressure
           xx=prs
        case(2)
           stop "eh 2"

        case(3)
           stop "eh 3"

        case(4)
           xx=eio

        case(5)
           stop "ehbah"

        case(6)
           xx=soundsqr
 
        case(7)
           xx=soundsqr*ri/prs

        case(8)
           xx=xx

        case(9)
           stop "eosflag=9 not implemented for hybrid eos"

        case default
           stop "eosflag not implemented for hybrid eos"

        end select
        
        return
     endif

  if(eoskey.eq.1) then

     select case (eosflag)

        case(1)
           !pressure
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=prs*press_gf
           tio=p_th
!           write(*,"(1P10E15.6)") t
        case(2)
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=dpdrho*press_gf/rho_gf

        case(3)
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=dpde*press_gf/eps_gf

        case(4)
           keytemp=1
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=e*eps_gf
           eio = xx

        case(5)
           stop "eosflag=5 not implemented for hybrid eos"

        case(6)
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           if (GR) then
              xx=soundsqr/(clite**2*(1.0d0+(prs*press_gf)/ri+eio))
           else
              xx=soundsqr/clite**2
           endif

        case(7)
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=soundsqr*r/prs

        case(8)
           stop "eosflag=8 not implemented for hybrid eos"

        case(9)
           stop "eosflag=9 not implemented for hybrid eos"

        case default
           stop "eosflag not implemented for hybrid eos"

     end select

  else if(eoskey.eq.2) then

     select case (eosflag)

        case(1)
           !pressure
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=prs*press_gf

        case(2)
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=dpdrho*press_gf/rho_gf

        case(3)
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=dpde*press_gf/eps_gf

        case(4)
           keytemp=1
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=e*eps_gf
           eio = xx

        case(5)
           stop "eosflag=5 not implemented for polytropic eos"

        case(6)
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           if (GR) then
              xx=soundsqr/(clite**2*(1.0d0+(prs*press_gf)/ri+eio))
           else
              xx=soundsqr/clite**2
           endif

        case(7)
           stop "eosflag=7 not implemented for polytropic eos"

        case(8)
           stop "eosflag=8 not implemented for polytropic eos"

        case(9)
           stop "eosflag=9 not implemented for polytropic eos"

        case default
           stop "eosflag not implemented for polytropic eos"

     end select

  else if(eoskey.eq.4) then

     select case (eosflag)

        case(1)
           !pressure
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=prs*press_gf

        case(2)
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=dpdrho*press_gf/rho_gf


        case(3)
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=dpde*press_gf/eps_gf

        case(4)
           keytemp=1
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=e*eps_gf
           eio = xx

        case(5)
           stop "eosflag=5 not implemented for ideal fluid"

        case(6)
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           if (GR) then
              xx=soundsqr/(clite**2*(1.0d0+(prs*press_gf)/ri+eio))
           else
              xx=soundsqr/clite**2
           endif

        case(7)
           stop "eosflag=7 not implemented for ideal fluid"

        case(8)
           stop "eosflag=8 not implemented for ideal fluid"

        case(9)
           stop "eosflag=9 not implemented for ideal fluid"

        case default
           stop "eosflag not implemented for ideal fluid"

     end select

  else if(eoskey.eq.3) then
#if HAVE_NUC_EOS
     select case (eosflag)
        case(1)
           !pressure 
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           xx=prs*press_gf
           if(keytemp.eq.1.or.keytemp.eq.2) then
              eio = e*eps_gf
           endif

        case(2)
           !dpdrhoe
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           xx=dpdrho*press_gf/rho_gf
           if(keytemp.eq.1.or.keytemp.eq.2) then
              eio = e*eps_gf
           endif

        case(3)
           !dpde
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           xx=dpde*press_gf/eps_gf
           if(keytemp.eq.1.or.keytemp.eq.2) then
              eio = e*eps_gf
           endif

        case(4)
           !eps
           keytemp=1
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           xx=e*eps_gf
           eio = xx

        case(5)
           stop "eosflag=5 not implemented for eoskey = 3"

        case(6)
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           if(keytemp.eq.1.or.keytemp.eq.2) then
              eio = e*eps_gf
           endif
           if (GR) then
              xx=soundsqr/(clite**2*(1.0d0+(prs*press_gf)/ri+eio))
           else
              xx=soundsqr/clite**2
           endif

        case(7)
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           xx = soundsqr*r/prs

        case(8)
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           if(keytemp.eq.1.or.keytemp.eq.2) then
              eio = e*eps_gf
           endif
           xx = xent

        case(9)
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           if (keytemp.eq.1.or.keytemp.eq.2) then
              eio = e*eps_gf
           endif
           xx = xmunu

        case default
           stop "eosflag not implemented for eoskey = 3"

        end select

        if(keyerr.eq.667) then
           if(i.ne.n1) then
              if(atmo(i+1).ne.0) then
                 atmo(i) = 1
                 keyerr = 0
                 call atmos_eos(i,ri,prs,eio,soundsqr)
                 select case (eosflag)
                 case(1)
                    !pressure
                    xx=prs
                 case(2)
                    stop "eh 4"
                    
                 case(3)
                    stop "eh 5"

                 case(4)
                    xx=eio
                    
                 case(5)
                    stop "ehbah"
                    
                 case(6)
                    xx=soundsqr
                    
                 case(7)
                    xx=soundsqr*ri/prs
                    
                 case(8)
                    xx=xx
                    
                 case(9)
                    stop "eosflag=9 not implemented for hybrid eos"
                    
                 case default
                    stop "eosflag not implemented for hybrid eos"
                    
                 end select
              endif
           endif
        endif

#endif

  else
     write(6,*) "eoskey: ", eoskey
     call flush(6)
     stop "eos choice not implemented, sorry..."

  endif

end subroutine eos

module hybrid_eos_module

  real*8 gamma1
  real*8 gamma2
  real*8 gammath
  real*8 K1,K2,E1,E2,E3
  real*8 rhonuc

end module hybrid_eos_module

module pwp_module
  implicit none
  integer, parameter :: npwp = 41
  real*8, parameter :: pwp_K(41) = (/ &
     6.1187577190167852d+12, &
     9.5501074299208838d+14, &
     4.7721408896840747d+22, &
     3.6573561788633219d+13, &
     3.4228609521189687d-04, &
     7.5169526027254232d+00, &
     8.7822319172585883d-61, &
     3.8807556960320082d-43, &
     2.7500152550046680d+34, &
     4.3637802033757011d-98, &
     3.4355024086293046d-93, &
     6.0411494828146495d-91, &
     3.2285009984275828d-88, &
     6.0518084864430667d-85, &
     4.4647244880479175d-81, &
     1.3811123638710344d-76, &
     1.7659778902332391d-71, &
     8.3253220844144581d-66, &
     1.1537582355309171d-59, &
     3.3675125934911839d-53, &
     1.3674472192504304d-46, &
     4.9154478684752059d-40, &
     1.0102461348657202d-33, &
     8.1727991488188082d-28, &
     1.9758969750446468d-22, &
     1.2119449801603048d-17, &
     1.7806724023668426d-13, &
     6.4590191422668059d-10, &
     6.3471752409492706d-07, &
     1.9260573505894103d-04, &
     2.0927219884958683d-02, &
     9.4564588128639626d-01, &
     2.0469945314776766d+01, &
     2.4110279502052703d+02, &
     1.7266416141150153d+03, &
     8.2625551106359544d+03, &
     2.8586306453008368d+04, &
     7.6294742563357737d+04, &
     1.6559619401600913d+05, &
     3.0502418827631086d+05, &
     3.0502418827631086d+05 /)
  real*8, parameter :: pwp_G(41) = (/ &
     1.5840000000000001d+00, &
     1.2869999999999999d+00, &
     6.2200000000000000d-01, &
     1.3560000000000001d+00, &
     2.5699999999999998d+00, &
     2.2697889258168567d+00, &
     6.4374483852130124d+00, &
     5.2415665266805638d+00, &
     4.8185899156652390d-02, &
     8.7973859713310070d+00, &
     8.4724499432229639d+00, &
     8.3234602614035378d+00, &
     8.1424506248027715d+00, &
     7.9253016472360178d+00, &
     7.6687107911884684d+00, &
     7.3708767025005351d+00, &
     7.0322418109642806d+00, &
     6.6561394156628264d+00, &
     6.2491393517368410d+00, &
     5.8208950710029352d+00, &
     5.3833946164812962d+00, &
     4.9496948201159823d+00, &
     4.5324024413164832d+00, &
     4.1422642294464334d+00, &
     3.7871848228614882d+00, &
     3.4718315796415538d+00, &
     3.1977945559025995d+00, &
     2.9641338406123645d+00, &
     2.7681026002870737d+00, &
     2.6058666869941365d+00, &
     2.4731096051982808d+00, &
     2.3654785466351282d+00, &
     2.2788743230265363d+00, &
     2.2096125195887684d+00, &
     2.1544905367601284d+00, &
     2.1107926489822186d+00, &
     2.0762585485055989d+00, &
     2.0490335493845184d+00, &
     2.0276123501688832d+00, &
     2.0107834982720152d+00, &
     2.0107834982720152d+00 /)
  real*8, parameter :: pwp_a(41) = (/ &
     0.0000000000000000d+00, &
     -2.2282553734253792d+17, &
     1.2081426138038909d+19, &
     6.8202175648473037d+18, &
     1.4648172077930039d+19, &
     1.2032852421530307d+19, &
     2.8686294330900783d+19, &
     2.0954829745689338d+19, &
     3.1201768525852869d+20, &
     1.6324984620610519d+20, &
     1.6249640014431725d+20, &
     1.6212185311448603d+20, &
     1.6163490917733889d+20, &
     1.6099961130715434d+20, &
     1.6016740088584405d+20, &
     1.5907230523066011d+20, &
     1.5762427862520375d+20, &
     1.5570014527496785d+20, &
     1.5313163177916419d+20, &
     1.4969026427251753d+20, &
     1.4506965433394784d+20, &
     1.3886716147123211d+20, &
     1.3056926015470846d+20, &
     1.1954794227063138d+20, &
     1.0507816315532331d+20, &
     8.6386678155336385d+19, &
     6.2738029349953430d+19, &
     3.3552421650188906d+19, &
     -1.4652023923345828d+18, &
     -4.2218651310002856d+19, &
     -8.8172103898080477d+19, &
     -1.3838122636129814d+20, &
     -1.9159962121034708d+20, &
     -2.4643202073718548d+20, &
     -3.0149421106561758d+20, &
     -3.5554468768541940d+20, &
     -4.0756876237143592d+20, &
     -4.5681298182551444d+20, &
     -5.0277956430170849d+20, &
     -5.4519524732403050d+20, &
     -5.4519524732403050d+20 /)
  real*8, parameter :: pwp_rj(40) = (/ &
     2.4266100950824142d+07, &
     3.7757219092541602d+11, &
     2.6242185433844443d+12, &
     1.0641430182243142d+14, &
     2.8973435877013212d+14, &
     4.1715867134425562d+14, &
     5.6893462243467100d+14, &
     6.2771657692207712d+14, &
     1.1592417091825145d+15, &
     1.1694274427643725d+15, &
     1.1722643780672420d+15, &
     1.1758861659928618d+15, &
     1.1805066484533650d+15, &
     1.1863959155021740d+15, &
     1.1938938949524770d+15, &
     1.2034265263584530d+15, &
     1.2155246171837282d+15, &
     1.2308452459379325d+15, &
     1.2501952245341620d+15, &
     1.2745556639911760d+15, &
     1.3051061567916175d+15, &
     1.3432466244407465d+15, &
     1.3906146945027532d+15, &
     1.4490968309253678d+15, &
     1.5208325479270795d+15, &
     1.6082128680675620d+15, &
     1.7138763562963235d+15, &
     1.8407078688882342d+15, &
     1.9918458584522610d+15, &
     2.1707033120013328d+15, &
     2.3810054416929380d+15, &
     2.6268448807301130d+15, &
     2.9127532228182860d+15, &
     3.2437868038052910d+15, &
     3.6256246932910815d+15, &
     4.0646776315635190d+15, &
     4.5682077081086560d+15, &
     5.1444596187651410d+15, &
     5.8028052007281910d+15, &
     6.5539035983828440d+15 /)
end module pwp_module

subroutine pwp_cold(rho, pco, ecold, dpco)
  use pwp_module
  implicit none
  real*8, intent(in)  :: rho
  real*8, intent(out) :: pco, ecold, dpco
  integer :: i, seg
  seg = 1
  do i = 1, npwp-1
     if (rho .ge. pwp_rj(i)) seg = i+1
  enddo
  pco   = pwp_K(seg) * rho**pwp_G(seg)
  ecold = pwp_a(seg) + pwp_K(seg) * rho**(pwp_G(seg)-1.d0) / (pwp_G(seg)-1.d0)
  dpco  = pwp_G(seg) * pco / rho
end subroutine pwp_cold

subroutine hybrid_eos(rho,enr,prs,pth,dpdrho_o,dpde_o,soundsqr,keytemp)
  ! replaced: cold part = our 15-segment piecewise polytrope; thermal part = (gammath-1)*rho*(eps-eps_cold)
  use hybrid_eos_module
  implicit none
  real*8 rho,enr,prs,soundsqr
  integer keytemp
  real*8 dpdrho_o,dpde_o,pth
  real*8 pcold,ecold,dpcold,eth,dpth_drho,dpth_denr,dp_drho,dpde
  call pwp_cold(rho, pcold, ecold, dpcold)
  if(keytemp.eq.1) then
     enr = ecold
     prs = pcold
     soundsqr = dpcold
     pth = 0.d0
     dpdrho_o = dpcold
     dpde_o = 0.d0
     return
  endif
  eth = rho*(enr - ecold)
  pth = (gammath - 1.d0)*eth
  dpth_drho = (gammath - 1.d0)*(enr - ecold)
  dpth_denr = (gammath - 1.d0)*rho
  if (pth.lt.0.0d0) then
     pth = 0.0d0; dpth_drho = 0.0d0; dpth_denr = 0.0d0
  endif
  prs = pcold + pth
  dpde = dpth_denr
  dp_drho = dpcold + dpth_drho
  soundsqr = dp_drho + dpde*prs/rho**2
  dpdrho_o = dp_drho
  dpde_o = dpde
end subroutine hybrid_eos



subroutine init_hybrid_eos

  use hybrid_eos_module

  implicit none

  E1 = K1/(gamma1-1.d0)
  E2 = (gamma1 - 1.d0)/(gamma2-1.d0)*E1*rhonuc**(gamma1-gamma2)
  K2 = (gamma2 - 1.d0)*E2
  E3 = (gamma2 - gamma1)/(gamma2-1.d0)*E1*rhonuc**(gamma1-1.d0)

end subroutine init_hybrid_eos

module poly_eos_module

  real*8 polygamma
  real*8 polyK

end module poly_eos_module

subroutine poly_eos(rho,enr,prs,dpdrho,dpde,soundsqr,keytemp)

  use poly_eos_module
  implicit none

! Input/Output
  real*8 rho,enr,prs,soundsqr
  integer keytemp
  real*8 dpdrho,dpde

! Local 
!
  prs = polyK * rho**polygamma
  soundsqr = polygamma*polyK*rho**(polygamma-1.d0)

  dpdrho = polyK * rho**(polygamma - 1) * polygamma
  dpde = 0.0d0

  if(keytemp.eq.1) then
!	energy wanted
     enr=prs/rho/(polygamma-1.d0)
     soundsqr=polygamma*prs/rho
  endif


end subroutine poly_eos

module ideal_eos_module

  real*8 idealgamma
  real*8 idealK1

end module ideal_eos_module

subroutine ideal_eos(rho,enr,prs,dpdrho,dpde,soundsqr,keytemp)

  use ideal_eos_module
  implicit none

! Input/Output
  real*8 rho,enr,prs,soundsqr
  real*8 dpdrho,dpde
  integer keytemp


  if(keytemp.eq.1) then
!	energy wanted
     prs=idealK1*rho**(idealgamma)
     enr=prs/rho/(idealgamma-1.d0)
     soundsqr=idealgamma*prs/rho

  endif
     
  prs = (idealgamma - 1.0d0) *rho*enr

  dpde = (idealgamma - 1.0d0 ) * rho
  dpdrho = (idealgamma - 1.0d0 ) * enr

  soundsqr= dpdrho+dpde*prs/rho**2


end subroutine ideal_eos

subroutine atmos_eos(i,xrho,xprs,xenr,xcs2)

  use atmos
  use GR1D_module, only: rho_gf,press_gf, &
       v1,v,W,atmo
  implicit none

  integer, intent(in) :: i
  real*8 :: xrho,xprs, xenr, xcs2
  real*8,parameter :: idealK1 =  1.2435d15 * (0.5d0**(4.d0/3.d0))
  real*8,parameter :: idealgamma = 1.66666666666d0

  xrho = atmo_rho
  xprs=idealK1*((atmo_rho/rho_gf)**(idealgamma))*press_gf
  xenr=xprs/atmo_rho/(idealgamma-1.d0)
  xcs2=idealgamma*xprs/atmo_rho
  v(i) = 0.0d0
  v1(i) = 0.0d0
  W(i) = 1.0d0
  atmo(i) = 1

end subroutine atmos_eos


!*************************************************************************
! Parametrised energy sink (stand-in for neutrino cooling), driven by environment variables:
!   COOL_TAU     e-folding time of the thermal energy [s]   (<=0 or unset: off)
!   COOL_RHO     only zones with rho > COOL_RHO [g/cm3]      (default 1e12)
!   COOL_TRIGGER start when max(rho) first exceeds this [g/cm3] (default 5e14)
!   COOL_ECAP    stop after this much energy at infinity is radiated [Msun c^2] (default huge)
! Thermal energy = eps - eps_cold(rho); it relaxes toward 0 as exp(-t/tau). Cold energy is never removed.
subroutine apply_cooling(dts)
  use GR1D_module
  implicit none
  real*8 dts
  logical, save :: init = .false., active = .false., finished = .false., enabled = .false.
  real*8, save :: tau = -1.d0, rhothr = 1.d12, trig = 5.d14, ecap = 1.d30, erad = 0.d0
  character*64 buf
  integer st, i, keyerr, keytemp
  real*8 pc, ecold, dp, eth, deps, dm, rcgs, fac, esum_th
  real*8 :: rv
  if (.not.init) then
     init = .true.
     call get_environment_variable('COOL_TAU', buf, status=st)
     if (st.eq.0) then
        read(buf,*) tau
        if (tau.gt.0.d0) enabled = .true.
     endif
     call get_environment_variable('COOL_RHO', buf, status=st)
     if (st.eq.0) read(buf,*) rhothr
     call get_environment_variable('COOL_TRIGGER', buf, status=st)
     if (st.eq.0) read(buf,*) trig
     call get_environment_variable('COOL_ECAP', buf, status=st)
     if (st.eq.0) read(buf,*) ecap
     if (enabled) then
        open(unit=777,file=trim(adjustl(outdir))//'/cooling_t.dat',status='unknown',form='formatted')
        write(*,*) 'COOLING ON: tau[s], rho_thr, trigger, ecap = ', tau, rhothr, trig, ecap
     endif
  endif
  if (.not.enabled .or. finished) return
  if (.not.active) then
     if (maxval(rho(ghosts1+1:n1-ghosts1))/rho_gf .gt. trig) active = .true.
     if (.not.active) return
  endif
  fac = 1.d0 - exp(-dts/(tau*time_gf))
  esum_th = 0.d0
  do i = ghosts1+1, n1-ghosts1
     rcgs = rho(i)/rho_gf
     if (rcgs .gt. rhothr) then
        call pwp_cold(rcgs, pc, ecold, dp)
        eth = eps(i) - ecold*eps_gf
        if (eth .gt. 0.d0) then
           deps = eth*fac
           dm = 4.d0*pi*x1(i)**2*X(i)*W(i)*rho(i)*(x1i(i+1)-x1i(i))
           erad = erad + alp(i)*deps*dm
           eps(i) = eps(i) - deps
           keytemp = 0; keyerr = 0
           call eos(i,rho(i),temp(i),ye(i),eps(i),press(i),keytemp,keyerr,1,eoskey,eos_rf_prec)
           call eos(i,rho(i),temp(i),ye(i),eps(i),cs2(i),keytemp,keyerr,6,eoskey,eos_rf_prec)
        endif
     endif
  enddo
  write(777,'(1P4E16.8)') time, erad, maxval(rho(ghosts1+1:n1-ghosts1))/rho_gf, dts/time_gf
  if (erad .gt. ecap) finished = .true.
  call prim2con
  call boundaries(0,0)
  if (GR) then
     call con2GR
     call GR_alp
     call GR_boundaries
  endif
end subroutine apply_cooling
