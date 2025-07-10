{$A+,B-,D+,E-,F+,G+,I+,L+,N+,O+,P+,Q+,R+,S+,T-,V+,X+,Y+}
{$M 65520,0,655360}
unit mpdq ;
interface

         uses mpvd,mnk,sbop;

         procedure cirtpvd(ro0,teta0,rot,tetat,krvz:extended;var vr,massa:extended;t,pia:extended; var ti:das);

         procedure preor(prd,abp,ai,omb,th,ev,omm:extended; nt:integer;t0:extended;et,ro,teta:dst);
         procedure tpvd(ro0,teta0,rot,tetat,krvz,vr,massa,t,pia:extended; var ti:das);
         procedure rtti (ti:das;t:extended; var ro,teta:extended);
         procedure rtti2(ti:das;t:extended; var x,y,ebc,cebc,sebc,xxc,yyc:extended);
         procedure rtti3(t,ev,th,prd,aa,bb,ff,gg:extended;var ro,teta:extended);
         procedure rtti4(ti:das;t:extended;var ro,teta:extended);


         {procedure vydres(ro0,teta0,rot,tetat,krvz,vr,massa,pia,t0:extended; et,ro,teta:dst;nt:integer;
         sgm:stb;bro0,bteta0,brot,btetat,br,bvr,bmassa,bpia:boolean;nro0,sro0,nteta0,steta0,nrot,srot,ntetat,stetat
         ,oevr,oevt:extended);}
         procedure maxnev(ro0,teta0,rot,tetat,krvz,vr,massa,pia,t0:extended; et,ro,teta:dst;nt:integer;
         var oevr,oevt:extended);
         {procedure nevotap(ro0,teta0,rot,tetat,krvz,vr,massa,pia,t0:extended; et,ro,teta:dst;nt:integer );}
         procedure vydres2(ro0,teta0,rot,tetat,krvz,vr,massa,pia,t0:extended; et,ro,teta:dst;nt:integer; oevr,oevt:extended );
         procedure gen0(ro0,teta0,rot,tetat,krvz,vr,massa,pia,t0:extended; et,ro,teta:dst;nt:integer; oevr,oevt:extended );
         procedure anres(ro0,teta0,rot,tetat,krvz,vr,massa,pia,t0:extended; et,ro,teta:dst;nt:integer; oevr,oevt,ev:extended );
         procedure vydres3(t0,pia:extended);
         procedure linprib(var xa,xt,xtt,ya,yt,ytt:extended;t1,nnt:integer;a:msls;cv:matrls;xf,sgm:stb;x,y,et:dst3);
         procedure linprib2(var xa,xt,xtt,ya,yt,ytt:extended;t1,nnt:integer;a:msls;cv:matrls;xf,sgm:stb;x,y,et:dst3;
                   var r,om,om0,xk,yk:extended;var pos:boolean;gbop:idst3;var ae6:ext6);
         procedure DrgEphPvd(iro0,iteta0,rot,tetat,r,vr,massa,t2:extended;
                             var diro0,diteta0,drot,dtetat,dr,dvr,dmassa:extended);
         procedure GetPrmEl(ev,th,prd:extended;var aa,bb,ff,gg,oev:extended);
         procedure GetElFrTi(var abp,ai,omb,omm:extended; aa,bb,ff,gg:extended);
         procedure GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,t2:extended;var iro0,iteta0,rot,tetat,r,vr,massa:extended);
         function GetUgol(abp1,ev1,th1,prd1,ai1,omb1,omm1,
                          abp2,ev2,th2,prd2,ai2,omm2,omb2:extended):extended;
         function GetUgol2(abp1,ev1,th1,prd1,ai1,omb1,omm1,
                          abp2,ev2,th2,prd2,ai2,omm2,omb2:extended):extended;
         procedure vydnev(tb:integer);






implementation


procedure linprib(var xa,xt,xtt,ya,yt,ytt:extended;t1,nnt:integer;a:msls;cv:matrls;xf,sgm:stb;x,y,et:dst3);
var rse:text;
xc,yc:extended;
km:integer;
begin
           km:=0;
           for i:=t1 to nnt do
           begin
                km:=km+1;
                a^[km,1]:=et^[i]*et^[i];
                a^[km,2]:=et^[i];
                a^[km,3]:=1;
                a^[km,4]:=x^[i];
           end;
           ptlst(a,xf,cv,km,3);
           sgmnk(a,cv,xf,sgm,3,km);
           xa:=xf[3];
           xt:=xf[2];
           xtt:=xf[1];
           km:=0;
           for i:=t1 to nnt do
           begin
                km:=km+1;
                a^[km,1]:=et^[i]*et^[i];
                a^[km,2]:=et^[i];
                a^[km,3]:=1;
                a^[km,4]:=y^[i];
           end;
           ptlst(a,xf,cv,km,3);
           sgmnk(a,cv,xf,sgm,3,km);
           ya:=xf[3];
           yt:=xf[2];
           ytt:=xf[1];
           {assign(rse,'linp.res');
           rewrite(rse);
           for i:=t1 to nnt do
           begin
                xc:=et[i]*et[i]*xtt+et[i]*xt+xa;
                yc:=et[i]*et[i]*ytt+et[i]*yt+ya;
                writeln(rse,et[i]:10:5,' ',x[i]:10:5,' ',y[i]:10:5);
                writeln(rse,et[i]:10:5,' ',xc:10:5,' ',yc:10:5);
           end;
           close(rse);}

end;
procedure linprib2(var xa,xt,xtt,ya,yt,ytt:extended;t1,nnt:integer;a:msls;cv:matrls;xf,sgm:stb;x,y,et:dst3;
                   var r,om,om0,xk,yk:extended;var pos:boolean;gbop:idst3;var ae6:ext6);
                       {mr,mom,mom0,mx0,my0}
var rse:text;
xc,yc,aa,bb,h,r4,mu,med,dxw,dyw,tom0,tom,tx0,ty0,ted,nevx,nevy:extended;
km,qt,l1,l2,j,q1,q2,ix:integer;
mj:longint;
kk,cc,tc1,tc2,tc3,x1,y1,y2,x2,xs,ys,kk2,cc2,a1,b1,c1,k1,c3,d1,xk1,xk2,yk1,yk2,r5
,r6,r7,r8,o1,o2,tt1,sg1,sg2,sg1q,sg2q,dlp,dlx,dlxt:extended;
{om0,om,xk,yk,r:extended;}
jc:integer;
xfx:stb;
btc1,btc2:boolean;

begin
           assign(rse,'linp.res');
           rewrite(rse);

           km:=0;
           for i:=t1 to nnt do
           begin
                if gbop^[i]=0 then
                begin
                     km:=km+1;
                     a^[km,1]:=et^[i];
                     a^[km,2]:=1;
                     a^[km,3]:=x^[i];
                end;
           end;
           ptlst(a,xf,cv,km,2);
           sgmnk(a,cv,xf,sgm,2,km);
           sg1:=sgm[3];
           xa:=xf[2];
           xt:=xf[1];
           xtt:=0;
           km:=0;
           q1:=32000;
           q2:=-32000;
           for i:=t1 to nnt do
           begin
                if gbop^[i]=0 then
                begin
                     km:=km+1;
                     a^[km,1]:=et^[i];
                     a^[km,2]:=1;
                     a^[km,3]:=y^[i];
                     if i<q1 then q1:=i;
                     if i>q2 then q2:=i;
                end;
           end;
           ptlst(a,xf,cv,km,2);
           sgmnk(a,cv,xf,sgm,2,km);
           sg2:=sgm[3];
           sg1:=sqrt(sg1*sg1+sg2*sg2);

           ya:=xf[2];
           yt:=xf[1];
           ytt:=0;

           km:=0;
           for i:=t1 to nnt do
           begin
                if gbop^[i]=0 then
                begin
                     km:=km+1;
                     a^[km,1]:=et^[i]*et^[i];
                     a^[km,2]:=et^[i];
                     a^[km,3]:=1;
                     a^[km,4]:=x^[i];
                end;
           end;
           ptlst(a,xfx,cv,km,3);
           sgmnk(a,cv,xfx,sgm,3,km);
           sg1q:=sgm[4];
           {writeln(sg1:10:5,' ',sgm[4]:10:5);}
           km:=0;
           for i:=t1 to nnt do
           begin
                if gbop^[i]=0 then
                begin
                     km:=km+1;
                     a^[km,1]:=et^[i]*et^[i];
                     a^[km,2]:=et^[i];
                     a^[km,3]:=1;
                     a^[km,4]:=y^[i];
                end;
           end;
           ptlst(a,xf,cv,km,3);
           sgmnk(a,cv,xf,sgm,3,km);
           {writeln(sg2:10:5,' ',sgm[4]:10:5);}
           sg2q:=sgm[4];
           sg1q:=sqrt(sg1q*sg1q+sg2q*sg2q);
           {writeln(sg1:10:5,' ',sg1q:10:5,' ',sg1/sg1q:10:5);}
           tc1:=et^[q1]+(et^[q2]-et^[q1])*0.25;
           tc2:=et^[q1]+(et^[q2]-et^[q1])*0.75;

           if (sg1q*2<sg1) {and false} then
           begin
                 xa:=xfx[3];
                 xt:=xfx[2];
                 xtt:=xfx[1];

                 ya:=xf[3];
                 yt:=xf[2];
                 ytt:=xf[1];
                 {writeln('ku-ku');
                 readln;}
                 aa:=et^[q1];
                 bb:=et^[q2];;
                 h:=(bb-aa)/100;
                 dlx:=0;
                 for ix:=0 to 100 do
                 begin
                      r4:=aa+ix*h;
                      xc:=r4*r4*xtt+r4*xt+xa;
                      yc:=r4*r4*ytt+r4*yt+ya;

                      if ix>0 then
                      begin
                           r5:=sqrt(sqr(xc-x1)+sqr(yc-y1));
                           dlx:=dlx+r5;
                      end;
                      x1:=xc;
                      y1:=yc;
                 end;
                 btc1:=false;
                 btc2:=false;
                 dlxt:=0;
                 for ix:=0 to 100 do
                 begin
                      r4:=aa+ix*h;
                      xc:=r4*r4*xtt+r4*xt+xa;
                      yc:=r4*r4*ytt+r4*yt+ya;

                      if ix>0 then
                      begin
                           r5:=sqrt(sqr(xc-x1)+sqr(yc-y1));
                           dlxt:=dlxt+r5;
                           if (dlxt>0.25*dlx) and not(btc1) then
                           begin
                                btc1:=true;
                                tc1:=r4;
                           end;
                           if (dlxt>0.75*dlx) and not(btc2) then
                           begin
                                btc2:=true;
                                tc2:=r4;
                           end;
                      end;
                      x1:=xc;
                      y1:=yc;
                 end;
           end;
           aa:=0;
           bb:=300;
           h:=0.1;
           qt:=trunc((bb-aa)/h);
           {mu:=sqrt(yt*yt+xt*xt);}

           {l1:=trunc((t1-1)+((nnt-t1+1)*0.25));
           if l1<t1 then l1:=t1;
           l2:=trunc((t1-1)+((nnt-t1+1)*0.75));
           if l2>nnt then l2:=nnt;}

           x1:=xtt*tc1*tc1+tc1*xt+xa;
           y1:=ytt*tc1*tc1+tc1*yt+ya;
           x2:=xtt*tc2*tc2+tc2*xt+xa;
           y2:=ytt*tc2*tc2+tc2*yt+ya;
           if x2<>x1 then
           begin
                kk:=(y2-y1)/(x2-x1);
           end else
           begin
               kk:=1;
           end;
           cc:=y1-kk*x1;
           {kk:=yt/xt;
           cc:=ya-yt*xa/xt;}


           ae6[1]:=x1;
           ae6[2]:=y1;
           ae6[3]:=tc1;
           ae6[4]:=x2;
           ae6[5]:=y2;
           ae6[6]:=tc2;
           tc3:=(tc2+tc1)/2;
           xs:=(x1+x2)/2;
           ys:=(y1+y2)/2;
           kk2:=-1/kk;
           cc2:=ys-kk2*xs;
           c3:=cc2-y1;
           FOR Jc:=0 To 100 do
           begin
                tt1:=jc*0.01;
           writeln(rse,0.0:10:5,' ',x1:10:5,y1:10:5);
           writeln(rse,0.0:10:5,' ',x1+tt1:10:5,y1:10:5);
           writeln(rse,0.0:10:5,' ',x1:10:5,y1+tt1:10:5);
           writeln(rse,0.0:10:5,' ',x1-tt1:10:5,y1:10:5);
           writeln(rse,0.0:10:5,' ',x1:10:5,y1-tt1:10:5);
           writeln(rse,0.0:10:5,' ',x2:10:5,y2:10:5);
           writeln(rse,0.0:10:5,' ',x2+tt1:10:5,y2:10:5);
           writeln(rse,0.0:10:5,' ',x2:10:5,y2+tt1:10:5);
           writeln(rse,0.0:10:5,' ',x2-tt1:10:5,y2:10:5);
           writeln(rse,0.0:10:5,' ',x2:10:5,y2-tt1:10:5);
           end;

           med:=1e39;
           pos:=true;
           for j:=0 to qt do
           begin
                r5:=sqrt(sqr(x2-x1)+sqr(y2-y1));
                r4:=abs(aa+j*h);
                a1:=1+kk2*kk2;
                b1:=-2*x1+2*kk2*c3;
                c1:=x1*x1+c3*c3-r4*r4;
                k1:=b1/2;
                d1:=(k1*k1-a1*c1);
                if d1>0 then
                begin
                     d1:=sqrt(d1);
                     xk1:=(-k1-d1)/a1;
                     xk2:=(-k1+d1)/a1;
                     yk1:=kk2*xk1+cc2;
                     yk2:=kk2*xk2+cc2;
                     r5:=sqrt(xk1*xk1+yk1*yk1);
                     r6:=sqrt(xk2*xk2+yk2*yk2);

                     r5:=sqrt(sqr(xk1-x1)+sqr(yk1-y1));
                     r6:=sqrt(sqr(xk1-x2)+sqr(yk1-y2));
                     r7:=sqrt(sqr(xk2-x1)+sqr(yk2-y1));
                     r8:=sqrt(sqr(xk2-x2)+sqr(yk2-y2));

                     dxw:=x2-xk1;
                     dyw:=y2-yk1;
                     o1:=arctan(dxw/dyw);
                     if dyw<0 then o1:=pi+o1;
                     if o1<0 then o1:=o1+2*pi;
                     o2:=o1;

                     dxw:=x1-xk1;
                     dyw:=y1-yk1;
                     o1:=arctan(dxw/dyw);
                     if dyw<0 then o1:=pi+o1;
                     if o1<0 then o1:=o1+2*pi;
                     tom:=(o2-o1)/(tc2-tc1);
                     tom0:=o1-tom*tc1;

                     {r5:=r4*sin(tom*tc1+tom0)+xk1;
                     r6:=r4*cos(tom*tc1+tom0)+yk1;
                     r7:=r4*sin(tom*tc2+tom0)+xk1;
                     r8:=r4*cos(tom*tc2+tom0)+yk1;}
                     ted:=0;
                     for i:=1 to nnt do
                     begin
                          if gbop^[i]=0 then
                          begin

                               xc:=r4*sin(tom*et^[i]+tom0)+xk1;
                               yc:=r4*cos(tom*et^[i]+tom0)+yk1;
                               nevx:=x^[i]-xc;
                               nevy:=y^[i]-yc;
                               ted:=ted+nevx*nevx+nevy*nevy;
                          end;
                     end;
                     ted:=sqrt(ted/(km*2));
                     if (ted<med) {and (r5<r6)} then
                     begin
                          med:=ted;
                          om0:=tom0;
                          om:=tom;
                          xk:=xk1;
                          yk:=yk1;
                          r:=r4;
                          pos:=r5<r6;
                          mj:=j;
                     end;
                     {if r4<20 then
                     begin
                          writeln(rse,0.0:10:5,' ',xk1:10:5,' ',yk1:10:5);
                     end;}
                     xk1:=xk2;
                     yk1:=yk2;

                     dxw:=x2-xk1;
                     dyw:=y2-yk1;
                     o1:=arctan(dxw/dyw);
                     if dyw<0 then o1:=pi+o1;
                     if o1<0 then o1:=o1+2*pi;
                     o2:=o1;

                     dxw:=x1-xk1;
                     dyw:=y1-yk1;
                     o1:=arctan(dxw/dyw);
                     if dyw<0 then o1:=pi+o1;
                     if o1<0 then o1:=o1+2*pi;
                     tom:=(o2-o1)/(tc2-tc1);
                     tom0:=o1-tom*tc1;

                     {r5:=r4*sin(tom*tc1+tom0)+xk1;
                     r6:=r4*cos(tom*tc1+tom0)+yk1;
                     r7:=r4*sin(tom*tc2+tom0)+xk1;
                     r8:=r4*cos(tom*tc2+tom0)+yk1;}
                     ted:=0;
                     for i:=1 to nnt do
                     begin
                          if gbop^[i]=0 then
                          begin
                               xc:=r4*sin(tom*et^[i]+tom0)+xk1;
                               yc:=r4*cos(tom*et^[i]+tom0)+yk1;
                               nevx:=x^[i]-xc;
                               nevy:=y^[i]-yc;
                               ted:=ted+nevx*nevx+nevy*nevy;
                          end;
                     end;
                     ted:=sqrt(ted/(km*2));
                     if (ted<med) {and (r5>=r6)} then
                     begin
                          med:=ted;
                          om0:=tom0;
                          om:=tom;
                          xk:=xk1;
                          yk:=yk1;
                          r:=r4;
                          pos:=r5>r6;
                          mj:=mj+10000;
                     end;
                     i:=0;
                end;
           end;
           ae6[7]:=r*sin(om*tc3+om0)+xk;
           ae6[8]:=r*cos(om*tc3+om0)+yk;
           ae6[9]:=tc3;

           ae6[7]:=tc3*tc3*xtt+tc3*xt+xa;;
           ae6[8]:=tc3*tc3*ytt+tc3*yt+ya;;
           ae6[9]:=tc3;

           i:=0;
           for i:=t1 to nnt do
           begin
                if gbop^[i]=0 then
                begin
                     xc:=r*sin(om*et^[i]+om0)+xk;
                     yc:=r*cos(om*et^[i]+om0)+yk;


                     writeln(rse,et^[i]:10:5,' ',x^[i]:10:5,' ',y^[i]:10:5);
                     writeln(rse,et^[i]:10:5,' ',xc:10:5,' ',yc:10:5);
                     xc:=et^[i]*et^[i]*xtt+et^[i]*xt+xa;
                     yc:=et^[i]*et^[i]*ytt+et^[i]*yt+ya;
                     writeln(rse,et^[i]:10:5,' ',xc:10:5,' ',yc:10:5);


                end;
           end;
           close(rse);
           ytt:=0;


end;




function eqkepler2(m,e:extended):extended;
const eps=1e-18;
var ep,eb,ep0:extended;
    i:integer;
begin
      ep:=m;
      i:=0;
      repeat
            i:=i+1;
            eb:=m+e*sin(ep);

            ep0:=ep;
            ep:=eb;
      until  (abs(ep0-eb)<eps) or (i>30000);
      eqkepler2:=eb;
      if i>30000 then
      begin
           writeln('Uravnenie keplera ne reshaetsa');
           readln;
           halt(0);
      end;
end;
function arsh(x:extended):extended;
begin
     arsh:=ln(x+sqrt(x*x+1));
end;
function hsin(x:extended):extended;
begin
     hsin:=0.5*(exp(x)-exp(-x));
end;

function hcos(x:extended):extended;
begin
     hcos:=0.5*(exp(x)+exp(-x));
end;
function eqkeplerhypold(m,e:extended):extended;
const eps=1e-12;
var ep,eb,ep0,de,fe:extended;
    i:integer;

begin
     i:=0;
      ep:=m;
      eb:=m;
      i:=0;
      repeat
            ep:=eb;
            eb:=arsh((eb+m)/e);
            i:=i+1;
      until  (abs(ep-eb)<eps) or (i>30000);
      eqkeplerhypold:=eb;
      if i>30000 then
      begin
           writeln('Uravnenie keplera ne reshaetsa3');
           readln;
           halt(0);
      end;

end;


function eqkeplerhyp(m,e:extended):extended;

const eps=1e-12;
      t30000=10000;
var ep,eb,ep0,de,fe,fh,fhh:extended;
    i:integer;

begin
     {if e>9 then
     begin
          i:=0;
     end;}

     i:=0;
      ep:=m;
      eb:=m;
      i:=0;
      repeat
            ep:=eb;
            fh:=e*hsin(eb)-eb-m;
            fhh:=e*hcos(eb)-1;
            eb:=eb-fh/fhh;
            i:=i+1;
      until  (abs(ep-eb)<eps) or (i>t30000);
      eqkeplerhyp:=eb;
      {if i>t30000 then
      begin
           writeln('Uravnenie keplera ne reshaetsa3');
           readln;
           halt(0);
      end;}

end;
function eqkepler(m,e:extended):extended;
const eps=1e-13;
      t30000=10000;
var ep,eb,ep0,de,fe,pi2:extended;
    i:integer;
begin
  pi2:=2*pi;
  while (m>pi) do m:=m-pi2;
  while (m<-pi) do m:=m+pi2;
  if m>10 then
  begin
       pi2:=0;
  end;
  if e>1 then
  begin
       {de:=eqkeplerhypold(m,e);}
       fe:=eqkeplerhyp(m,e);
       eqkepler:=fe;
  end else
  begin

      ep:=m;
      i:=0;
      repeat
            i:=i+1;
            fe:=ep-e*sin(ep)-m;
            de:=1-e*cos(ep);
            (*writeln(i,' ',fe,' ',de);
            if ((fe+6.11722142858814E+0017))<1 {and (de=1.12309171515122E-0002)} then
            begin
                  i:=0;
            end;*)

            eb:=ep-fe/de;
            while (eb>10*pi) do eb:=eb-pi2;
            while (eb<-10*pi) do eb:=eb+pi2;


            ep0:=ep;
            ep:=eb;
            i:=i+1;
            if abs(eb)>1e10 then i:=30001;
      until  (abs(ep0-eb)<eps) or (i>t30000);
      eqkepler:=eb;
      if i>t30000 then
      begin
           writeln('Uravnenie keplera ne reshaetsa2');
           {readln;
           halt(0);}

      ep:=m;
      i:=0;
      m:=m-pi2;
      repeat
            i:=i+1;
            fe:=ep-e*sin(ep)-m;
            de:=1-e*cos(ep);

            eb:=ep-fe/de;
            while (eb>100*pi) do eb:=eb-pi2;
            while (eb<-100*pi) do eb:=eb+pi2;


            ep0:=ep;
            ep:=eb;
            i:=i+1;
            if abs(eb)>1e10 then i:=30001;
      until  (abs(ep0-eb)<eps) or (i>t30000);
      end;


  end;
end;
procedure GetElFrTi(var abp,ai,omb,omm:extended; aa,bb,ff,gg:extended);
var rr,ss,cosai,ombMomm,dxw,dyw,ombPomm:extended;
begin
     rr:=0.5*(aa*aa+bb*bb+ff*ff+gg*gg);
     ss:=aa*gg-bb*ff;
     abp:=sqrt(rr+sqrt((rr+ss)*(rr-ss)));
     cosai:=(aa*gg-bb*ff)/(abp*abp);
     ai:=arccos(cosai);

     dxw:=bb+ff;
     dyw:=aa-gg;
     ombMomm:=arctan(dxw/dyw);
     if dyw<0 then ombMomm:=pi+ombMomm;
     if ombMomm<0 then ombMomm:=ombMomm+2*pi;

     dxw:=bb-ff;
     dyw:=aa+gg;
     ombPomm:=arctan(dxw/dyw);
     if dyw<0 then ombPomm:=pi+ombPomm;
     if ombPomm<0 then ombPomm:=ombPomm+2*pi;

     omb:=(ombPomm+ombMomm)*0.5;
     if omb<0 then omb:=omb+2*pi;
     if omb>2*pi then omb:=omb-2*pi;

     omm:=ombPomm-omb;
     if omm<0 then omm:=omm+2*pi;
     if omm>2*pi then omm:=omm-2*pi;

            omb:=pi/2-omb;
            if omb<0 then omb:=omb+2*pi;
            if omb>2*pi then omb:=omb-2*pi;


end;
procedure preor;
var dt0,rs3:text;
    comm,somm,comb,somb,ci,si,aa,bb,ff,gg,eb,xx,yy,xc,yc,rog,tetag,xo,yo
    ,seb,ceb,ad,hd,bd,etd:extended;
    i:integer;
begin
     assign(rs3,'traekt.res');
     append(rs3);

      th:=th-t0;
      assign(dt0,'preor.res');
      rewrite(dt0);
            ai:=ai+180;
            omm:=omm*pi/180;
            omb:=omb*pi/180;
            ai:=ai*pi/180;
            comm:=cos(omm);
            somm:=sin(omm);
            comb:=sin(omb);
            somb:=cos(omb);
            ci:=cos(ai);
            si:=sin(ai);
            aa:= comm*comb-somm*somb*ci;
            bb:= comm*somb+somm*comb*ci;
            ff:=-somm*comb-comm*somb*ci;
            gg:=-somm*somb+comm*comb*ci;
            for i:=1 to nt do
            begin
                 eb:=eqkepler((et[i]-th)*2*pi/prd,ev);

                 seb:=sin(eb);
                 ceb:=cos(eb);
                 xx:=ceb-ev;
                 yy:=sqrt(1-Ev*ev)*seb;
                 xc:=abp*(aa*xx+ff*yy) ;
                 yc:=abp*(bb*xx+gg*yy);
                 rog:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);

                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 tetag:=tetag*180/pi;
                 xo:=ro[i]*sin(teta[i]*pi/180);
                 yo:=ro[i]*cos(teta[i]*pi/180);

                 writeln(dt0,i:4,' ',(et[i]+t0):8:3,' ',ro[i]:8:3,' ',teta[i]:8:3,' ',(ro[i]-rog):8:3,' ',(teta[i]-tetag):8:3);

            end;
            close(dt0);
            ad:=et[1];
            bd:=et[nt];
            hd:=(bd-ad)/1000;

            for i:=0 to 1000 do
            begin
                 etd:=ad+hd*i;
                 eb:=eqkepler((etd-th)*2*pi/prd,ev);
                 seb:=sin(eb);
                 ceb:=cos(eb);
                 xx:=ceb-ev;
                 yy:=sqrt(1-Ev*ev)*seb;
                 xc:=abp*(aa*xx+ff*yy) ;
                 yc:=abp*(bb*xx+gg*yy);
                 rog:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);

                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 tetag:=tetag*180/pi;
                 writeln(rs3,etd+t0:8:3,' ',xc:20:10,' ',yc:20:10,' ',rog:20:10,' ',tetag:20:10);

            end;
            close(rs3);

end;

function por(x,p:extended):extended;
begin
    {writeln('x=',x:10:5,' p=',p:10:5);}
     if x=0 then
     begin
          por:=0;
     end else
     begin
          por:=exp(p*ln(x));
     end;

end;

procedure elor(var x,y,z,xt,yt,zt,akapa,a,e,omm,omb,ai,th,e0:mtreal);
var pf,r,v,rrt,se0,ce0:mtreal;
    px,py,pz,qx,qy,qz,ak,al,am,an, si,somm,somb,comm,comb :mtreal;
begin
      { Procedure для вычисления злементов орбиты по положению и скорости };
      {  Алгоритм взят из 'Введения в теор. астрономию' Субботина стр. 131 };
      {  x,y,z радиус вектор тела                 };
      {  xt,yt,zt вектор скорости                 };
      {  akapa=sqrt(m*g) для относительного движения    };
      {  a большая полуось                        };
      {  e экцентриситет                          };
      {  omm аргумент перицентра                  };
      {  omb долгота восходящего узла             };
      {  ai наклон орбиты                         };
      { readln(x,y,z,xt,yt,zt,amt,g  );}
      a:=0;e:=0;omm:=0;omb:=0;ai:=0;
      r:=sqrt(x*x+y*y+z*z);
      rrt:=x*xt+y*yt+z*zt;
      v:=xt*xt+yt*yt+zt*zt;
      a:=(2/r)-v/(akapa*akapa);

      a:=1/a;
     { writeln('a=', a ,a/r); }
      pf:=(v*r/(akapa*akapa))-1;
      e:= sqrt((rrt*rrt/(akapa*akapa*a))+(pf*pf));
                                                   if (e<>0.0) and (e<1) then
                                                   begin
      se0:=rrt/(akapa*sqrt(a)*e);
      ce0 :=pf/e;
      ak:=ce0/r;
      al:=sqrt(a)*se0/akapa;

      px:=ak*x-al*xt;
      py:=ak*y-al*yt;
      pz:=ak*z-al*zt;

      am:=se0/(r*sqrt(1-e*e));
      an:=sqrt(a/(1-e*e))*(ce0-e)/akapa ;
      qx:=am*x+an*xt;
      qy:=am*y+an*yt;
      qz:=am*z+an*zt;
      si:=sqrt(pz*pz+qz*qz);
      {writeln(otl2,se0:15:10,' ',ce0:15:10,' ',ak:15:10,' ',al:15:10,' ',
                am:15:10,' ',an:15:10,' ',si:15:10);
      flush(otl2);}

                                                         if (si<>0.0)  then
                                                         begin
      somm:=pz/si;
      comm:=qz/si;
      omm:=arctan(somm/sqrt(1-somm*somm));
      if comm<0 then  omm :=pi-omm;
      if omm<0 then omm:=2*pi+omm;

      somb:=py*comm-qy*somm;
      comb:=px*comm-qx*somm;
      omb:=arctan(somb/sqrt(1-somb*somb));
      if comb<0 then omb:=pi-omb;
      if omb<0 then omb:=2*pi+omb;
      {writeln(si);}
      if (abs(si-1))<0.0000001 then ai:=pi/2 else
      ai:=arctan(si/sqrt(1-si*si));
      if (x*yt-y*xt)<0 then ai:=pi-ai;

                                                         end;
                                                   end
                                                   else
                                                   begin
                                                         se0:=0;
                                                         ce0:=1;
                                                         ai:=0.5;


                                                         omb:=0.5;
                                                         omm:=0.5;
                                                   end;
      e0:=arcsin(se0);
      if ce0<0 then e0:=pi-e0;
      if e0<0 then e0:=2*pi+e0;
      a:=abs(a);
     th:=-(e0-e*se0)/(akapa*por(abs(a),-3/2));

end;

procedure elor3(var x,y,z,xt,yt,zt,akapa,a,e,omm,omb,ai,th,e0:mtreal);
var pf,r,v,rrt,se0,ce0:mtreal;
    px,py,pz,qx,qy,qz,ak,al,am,an, si,somm,somb,comm,comb,an2 :mtreal;
begin
      { Procedure для вычисления злементов орбиты по положению и скорости };
      {  Алгоритм взят из 'Введения в теор. астрономию' Субботина стр. 131 };
      {  x,y,z радиус вектор тела                 };
      {  xt,yt,zt вектор скорости                 };
      {  akapa=sqrt(m*g) для относительного движения    };
      {  a большая полуось                        };
      {  e экцентриситет                          };
      {  omm аргумент перицентра                  };
      {  omb долгота восходящего узла             };
      {  ai наклон орбиты                         };
      { readln(x,y,z,xt,yt,zt,amt,g  );}
      a:=0;e:=0;omm:=0;omb:=0;ai:=0;
      r:=sqrt(x*x+y*y+z*z);
      rrt:=x*xt+y*yt+z*zt;
      v:=xt*xt+yt*yt+zt*zt;
      if akapa<>0 then
      begin
            a:=(2/r)-v/(akapa*akapa);
            a:=1/a;
      end else
      begin
           a:=0;
      end;
     { writeln('a=', a ,a/r); }
      pf:=(v*r/(akapa*akapa))-1;
      e:= sqrt((rrt*rrt/(akapa*akapa*a))+(pf*pf));
      {e:=0;}


      if ((bcir) and (e<1e-6)) or (e=0) then
      begin
            se0:=0;
            ce0:=1;
            e:=0;

      end else
      begin
            se0:=rrt/(akapa*sqrt(abs(a))*e);
            ce0 :=pf/e;

      end;
      ak:=ce0/r;
      al:=sqrt(abs(a))*se0/akapa;




      px:=ak*x-al*xt;
      py:=ak*y-al*yt;
      pz:=ak*z-al*zt;

      am:=se0/(r*sqrt(abs(1-e*e)));
      an:= sqrt(a/(1-e*e))*(ce0-e)/akapa;
      if a<0 then an:=-an;
      {an:=abs(an);}


      qx:=am*x+an*xt;
      qy:=am*y+an*yt;
      qz:=am*z+an*zt;
      si:=sqrt(pz*pz+qz*qz);

      {writeln(otl2,se0:15:10,' ',ce0:15:10,' ',ak:15:10,' ',al:15:10,' ',
                am:15:10,' ',an:15:10,' ',si:15:10);
      flush(otl2);}

                                                         if (si<>0.0)  then
                                                         begin
      somm:=pz/si;
      comm:=qz/si;
      omm:=arctan(somm/sqrt(1-somm*somm));
      if comm<0 then  omm :=pi-omm;
      if omm<0 then omm:=2*pi+omm;

      somb:=py*comm-qy*somm;
      comb:=px*comm-qx*somm;
      omb:=arctan(somb/sqrt(1-somb*somb));
      if comb<0 then omb:=pi-omb;
      if omb<0 then omb:=2*pi+omb;
      {writeln(si);}
      if (abs(si-1))<0.0000001 then ai:=pi/2 else
      ai:=arctan(si/sqrt(1-si*si));
      if (x*yt-y*xt)<0 then ai:=pi-ai;

                                                         end;
      if e<1 then
      begin
           e0:=arcsin(se0);
           if ce0<0 then e0:=pi-e0;
           if e0<0 then e0:=2*pi+e0;
      end else
      begin
           e0:=arsh(se0);
      end;
      {a:=abs(a);}
      th:=-(e0-e*se0)/(akapa*por(abs(a),-3/2));
      if e>=1 then th:=-th;


end;
procedure elor2(var x,y,z,xt,yt,zt,akapa,a,e,omm,omb,ai,th,e0:mtreal);
begin
     {elor (x,y,z,xt,yt,zt,akapa,a,e,omm,omb,ai,th,e0);}
     elor3(x,y,z,xt,yt,zt,akapa,a,e,omm,omb,ai,th,e0);
     omb:=pi/2-omb;
     while omb>2*pi do omb:=omb-2*pi;
     while omb<0    do omb:=omb+2*pi;
end;

function ctg(x:extended):extended;
begin
     ctg:=cos(x)/sin(x);
end;
procedure rttipar(ti:das;t:extended; var ro,teta:extended);
var eb,seb,ceb,tetag,xx,yy,xc,yc,abp,ev,th,prd,aa,bb,ff,gg,mm,q,b,tmassa
,b2,b3,b4,leps,s02,s0,seps,dze,sinsg,sg,gog,vdze,udze:extended;
itr:integer;
begin
                 abp:=ti[1];
                 ev:=ti[2];
                 th:=ti[3];
                 prd:=ti[4];
                 aa:=ti[5];
                 bb:=ti[6];
                 ff:=ti[7];
                 gg:=ti[8];
                 tmassa:=abp*abp*abp/(ti[4]*ti[4]);
                 eb:=0;
                 q:=abp*(1-ev);
                 b:=(t-th)/(q*q*sqrt(q));
                 kp:=(4*pi*pi)*tmassa;
                 b2:=sqrt(kp/2)*b;
                 s0:=0;
                 s02:=0;
                 itr:=0;
                 leps:=1e-16;

                 repeat
                       s0:=s02;
                       b3:=s0-b2+s0*s0*s0/3;
                       b4:=1+s0*s0*2/3;
                       s02:=s0-b3/b4;
                       itr:=itr+1;
                 until (abs(s0-s02)<leps) or (itr>30000);

                        seps:=0.5*(1-ev);
                        dze:=seps*s02*s02;
                        sinsg:=sqrt(dze);
                        sg:=arcsin(sinsg);
                        gog:=(2*sg-sin(2*sg))/(4*sinsg*sinsg*sinsg);
                        vdze:=gog*sqrt(2/kp);
                        udze:=sqrt(1-dze)*sqrt(2/kp);
                        eb:=0;



end;
procedure GetPrmElOsn(ev,th,prd:extended;var aa,bb,ff,gg,oev:extended);
var i,km:integer;
mm,eb,xx,yy:extended;
xf:stb;

begin
     km:=0;
     for i:=1 to nt do
     begin
           if gbop^[i]=0 then
           begin
                 mm:=(et^[i]-th)*2*pi/prd;
                 eb:=eqkepler(mm,ev);
                 if ev<1 then
                 begin
                      seb:=sin(eb);
                      ceb:=cos(eb);
                 end else
                 begin
                      seb:=hsin(eb);
                      ceb:=hcos(eb);
                 end;
                 if ev<1 then
                 begin
                      xx:=ceb-ev;
                 end else
                 begin
                      xx:=-ceb+ev;
                 end;
                 yy:=sqrt(abs(1-Ev*ev))*seb;
                 km:=km+1;
                 a^[km,1]:=xx*aves^[i];
                 a^[km,2]:=yy*aves^[i];;
                 a^[km,3]:=X^[i]*aves^[i];;

                 at^[km,1]:=xx*aves^[i];;
                 at^[km,2]:=yy*aves^[i];;
                 at^[km,3]:=y^[i]*aves^[i];;
           end;

     end;
     ptlst(a,xf,cv,km,2);
     sgmnk(a,cv,xf,sgm,2,km);
     oev:=sgm[3];
     aa:=xf[1];
     ff:=xf[2];

     ptlst(at,xf,cv,km,2);
     sgmnk(at,cv,xf,sgm,2,km);
     oev:=sqrt(oev*oev+sgm[3]*sgm[3])*0.5;
     bb:=xf[1];
     gg:=xf[2];

end;
procedure GetPrmElFixoi(ev,th,prd:extended;var aa,bb,ff,gg,oev:extended);
var i,km:integer;
mm,eb,xx,yy,c,s,cc,sc,comb,somb,ci,ll,kk:extended;
xf:stb;
abp2,ai2,omb2,omm2:extended;

begin
     comb:=sin(omb);
     somb:=cos(omb);
     ci:=cos(ai);
     c:=comb;
     s:=somb;
     cc:=comb*ci;
     sc:=somb*ci;

     km:=0;
     for i:=1 to nt do
     begin
           if gbop^[i]=0 then
           begin
                 mm:=(et^[i]-th)*2*pi/prd;
                 eb:=eqkepler(mm,ev);
                 if ev<1 then
                 begin
                      seb:=sin(eb);
                      ceb:=cos(eb);
                 end else
                 begin
                      seb:=hsin(eb);
                      ceb:=hcos(eb);
                 end;
                 if ev<1 then
                 begin
                      xx:=ceb-ev;
                 end else
                 begin
                      xx:=-ceb+ev;
                 end;
                 yy:=sqrt(abs(1-Ev*ev))*seb;
                 km:=km+1;
                 al^[km]^[1]:= xx*c -yy*sc;
                 al^[km]^[2]:=-xx*sc-yy*c;
                 al^[km]^[3]:=X^[i];

                 km:=km+1;
                 al^[km]^[1]:=xx*s +yy*cc;
                 al^[km]^[2]:=xx*cc-yy*s;
                 al^[km]^[3]:=y^[i];
           end;

     end;
     ptlst4(al,xf,cv,km,2);
     sgmnk2(al,cv,xf,sgm,2,km);
     oev:=sgm[3];
     ll:=xf[1];
     kk:=xf[2];
     {
            aa:= comm*comb-somm*somb*ci;
            bb:= comm*somb+somm*comb*ci;
            ff:=-somm*comb-comm*somb*ci;
            gg:=-somm*somb+comm*comb*ci;}

     aa:= ll*c-kk*sc;
     bb:= ll*s+kk*cc;
     ff:=-kk*c-ll*sc;
     gg:=-kk*s+ll*cc;
     GetElFrTi(abp2,ai2,omb2,omm2,aa,bb,ff,gg);
     i:=0;

end;

(*procedure GetPrmElFixoi(ev,th,prd:extended;var aa,bb,ff,gg,oev:extended);
var i,km:integer;
mm,eb,xx,yy,c,s,cc,sc,comb,somb,ci,ll,kk:extended;
xf:stb;
abp2,ai2,omb2,omm2:extended;

begin
     comb:=sin(omb);
     somb:=cos(omb);
     ci:=cos(ai);
     c:=comb;
     s:=somb;
     cc:=comb*ci;
     sc:=somb*ci;

     km:=0;
     for i:=1 to nt do
     begin
           if gbop^[i]=0 then
           begin
                 mm:=(et^[i]-th)*2*pi/prd;
                 eb:=eqkepler(mm,ev);
                 if ev<1 then
                 begin
                      seb:=sin(eb);
                      ceb:=cos(eb);
                 end else
                 begin
                      seb:=hsin(eb);
                      ceb:=hcos(eb);
                 end;
                 if ev<1 then
                 begin
                      xx:=ceb-ev;
                 end else
                 begin
                      xx:=-ceb+ev;
                 end;
                 yy:=sqrt(abs(1-Ev*ev))*seb;
                 km:=km+1;
                 a^[km,1]:= xx*c -yy*sc;
                 a^[km,2]:=-xx*sc-yy*c;
                 a^[km,3]:=X^[i];

                 km:=km+1;
                 a^[km,1]:=xx*s +yy*cc;
                 a^[km,2]:=xx*cc-yy*s;
                 a^[km,3]:=y^[i];
           end;

     end;
     ptlst(a,xf,cv,km,2);
     sgmnk(a,cv,xf,sgm,2,km);
     oev:=sgm[3];
     ll:=xf[1];
     kk:=xf[2];
     {
            aa:= comm*comb-somm*somb*ci;
            bb:= comm*somb+somm*comb*ci;
            ff:=-somm*comb-comm*somb*ci;
            gg:=-somm*somb+comm*comb*ci;}

     aa:= ll*c-kk*sc;
     bb:= ll*s+kk*cc;
     ff:=-kk*c-ll*sc;
     gg:=-kk*s+ll*cc;
     GetElFrTi(abp2,ai2,omb2,omm2,aa,bb,ff,gg);
     i:=0;

end;*)

procedure GetPrmEl;
begin
     if fixoi then
     begin
            GetPrmElFixoi(ev,th,prd,aa,bb,ff,gg,oev);

     end else
     begin
            GetPrmElOsn(ev,th,prd,aa,bb,ff,gg,oev);
     end;



end;

procedure rtti(ti:das;t:extended; var ro,teta:extended);

var eb,seb,ceb,tetag,xx,yy,xc,yc,abp,ev,th,prd,aa,bb,ff,gg,mm:extended;
eb2:extended;

begin
                 {t:=-193;}
                 abp:=ti[1];
                 ev:=ti[2];
                 th:=ti[3];
                 prd:=ti[4];
                 aa:=ti[5];
                 bb:=ti[6];
                 ff:=ti[7];
                 gg:=ti[8];
                 {if (0.99<=ev) and (ev<=1) then
                 begin
                      rttipar(ti,t,ro,teta);
                 end;}
                 mm:=(t-th)*2*pi/prd;
                 eb:=eqkepler(mm,ev);
                 if ev<1 then
                 begin
                      seb:=sin(eb);
                      ceb:=cos(eb);
                 end else
                 begin
                      seb:=hsin(eb);
                      ceb:=hcos(eb);
                 end;
                 if ev<1 then
                 begin
                      xx:=ceb-ev;
                 end else
                 begin
                      xx:=-ceb+ev;
                 end;
                 yy:=sqrt(abs(1-Ev*ev))*seb;
                 xc:=abs(abp)*(aa*xx+ff*yy) ;
                 yc:=abs(abp)*(bb*xx+gg*yy);
                 ro:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);
                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 teta:=tetag*180/pi;

                 {writeln(otl2,eb:15:10,' ',seb:15:10,' ',ceb:15:10,' ',xx:15:10,' ',
                 yy:15:10,' ',xc:15:10,' ',yc:15:10);
                 flush(otl2);
                 xc:=0;}


end;
procedure rtti4;
var
ev,th,prd,aa,bb,ff,gg:extended;
begin
                 ev:=ti[2];
                 th:=ti[3];
                 prd:=ti[4];
                 aa:=ti[5];
                 bb:=ti[6];
                 ff:=ti[7];
                 gg:=ti[8];

     rtti3(t,ev,th,prd,aa,bb,ff,gg,ro,teta);
end;

procedure rtti3;

var eb,seb,ceb,tetag,xx,yy,xc,yc,mm:extended;
eb2:extended;

begin
                 mm:=(t-th)*2*pi/prd;
                 eb:=eqkepler(mm,ev);
                 if ev<1 then
                 begin
                      seb:=sin(eb);
                      ceb:=cos(eb);
                 end else
                 begin
                      seb:=hsin(eb);
                      ceb:=hcos(eb);
                 end;
                 if ev<1 then
                 begin
                      xx:=ceb-ev;
                 end else
                 begin
                      xx:=-ceb+ev;
                 end;
                 yy:=sqrt(abs(1-Ev*ev))*seb;
                 xc:=(aa*xx+ff*yy) ;
                 yc:=(bb*xx+gg*yy);
                 ro:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);
                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 teta:=tetag*180/pi;

end;

procedure rtti2;
var eb,seb,ceb,tetag,xx,yy,xc,yc,abp,ev,th,prd,aa,bb,ff,gg,mm:extended;

begin
                 abp:=ti[1];
                 ev:=ti[2];
                 th:=ti[3];
                 prd:=ti[4];
                 aa:=ti[5];
                 bb:=ti[6];
                 ff:=ti[7];
                 gg:=ti[8];
                 mm:=(t-th)*2*pi/prd;
                 eb:=eqkepler(mm,ev);
                 if ev<1 then
                 begin
                      seb:=sin(eb);
                      ceb:=cos(eb);
                 end else
                 begin
                      seb:=hsin(eb);
                      ceb:=hcos(eb);
                 end;
                 if ev<1 then
                 begin
                      xx:=ceb-ev;
                 end else
                 begin
                      xx:=-ceb+ev;
                 end;
                 yy:=sqrt(abs(1-Ev*ev))*seb;
                 xc:=abs(abp)*(aa*xx+ff*yy) ;
                 yc:=abs(abp)*(bb*xx+gg*yy);
                 x:=xc;
                 y:=yc;
                 ebc:=eb;
                 cebc:=ceb;
                 sebc:=seb;
                 xxc:=xx;
                 yyc:=yy;

end;

procedure tpvd;

var xs,ys,vx,vy,dxw,dyw,psi,mu,kp,xa,ya,rt,za,xt,yt,zt,omm,omb,ai,e0,
comm,somm,comb,somb,ci,si,eb,sep,ceb,xx,yy,xc,yc,seb,tetag,rog,ev0,ro:extended;
 abp,ev,th,prd,aa,bb,ff,gg,omb2,zaf:extended;
begin
                {writeln(dt0,' krivizna=',krvz:10:5);
                writeln(dt0,'teta0=',(teta0):8:3);
                writeln(dt0,'  ro0=',ro0:8:3);   }
                xs:=ro0*sin(teta0*pi/180);
                ys:=ro0*cos(teta0*pi/180);
                vx:=rot*sin(teta0*pi/180)+ro0*tetat*cos(teta0*pi/180)*pi/180;
                vy:=rot*cos(teta0*pi/180)-ro0*tetat*sin(teta0*pi/180)*pi/180;
                teta0:=teta0*pi/180;

                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;

                mu:=sqrt(vx*vx+vy*vy);
                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                mu:=sqrt(vx*vx+vy*vy);
                {writeln(dt0,'mu=',mu:10:5,' psi=',psi*180/pi:10:5,' deg');}
                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                kp:=(4*pi*pi)*massa;
                xa:=xs/pia;
                ya:=ys/pia;
                {if (kp*ro0*krvz*abs(sin(psi-teta0)))/(mu*mu)<0 then
                begin
                     rt:=por((kp*ro0*krvz*abs(sin(psi-teta0)))/(mu*mu),1/3);
                end;}
                rt:=por(abs((kp*ro0*krvz*abs(sin(psi-teta0)))/(mu*mu)),1/3);
                ro:=ro0/pia;
                {writeln(dt0,'r=',rt:8:2,' a.e.  ro=',ro:8:2,' a.e.');
                writeln(dt0,' massa=',massa:4:1,' pi=',pia:8:4,'"');}
                za:=abs(rt*rt-ro*ro);
                if za>0 then
                begin
                     za:=sqrt(za);
                end else
                begin
                     za:=0;
                end;
                if oprza then  za:=-za;
                {writeln(' za=',za:8:2,' a.e.');}
                xt:=mu*sin(psi)/pia;
                yt:=mu*cos(psi)/pia;
                zt:=vr*(3600.0*24.0*365.2422)/149600000.0;
                {zt:=vr;}
                kp:=sqrt(kp);

                {writeln(otl2,xa:15:10,' ',ya:15:10,' ',za:15:10,' ',xt:15:10,' ',
                yt:15:10,' ',zt:15:10,' ',kp:15:10);
                flush(otl2);}
                zaf:=-za;
                elor2(xa,ya,za,xt,yt,zt,kp,abp,ev,omm,omb,ai,th,e0);
                {th:=-121.4407430020;}

                {pvdv (xa,ya,za,xt,yt,zt,kp,abp,ev,ai,omb,omm,e0,th);}
                {writeln(otl2,abp:25:10,' ',ev:15:10,' ',omm:15:10,' ',omb:15:10,' ',
                ai:15:10,' ',th:15:10,' ',e0:15:10);
                flush(otl2);}


                ti[2]:=ev;


                {if ev>1 then exit;}
                prd:=kp/(sqrt(abs(abp*abp*abp))*2*pi);
                prd:=1/prd;
                abp:=abp*pia;


            comm:=cos(omm);
            somm:=sin(omm);
            comb:=sin(omb);
            somb:=cos(omb);
            ci:=cos(ai);
            si:=sin(ai);
            aa:= comm*comb-somm*somb*ci;
            bb:= comm*somb+somm*comb*ci;
            ff:=-somm*comb-comm*somb*ci;
            gg:=-somm*somb+comm*comb*ci;

                 {eb:=eqkepler((t-th)*2*pi/prd,ev);
                 seb:=sin(eb);
                 ceb:=cos(eb);
                 xx:=ceb-ev;
                 yy:=sqrt(abs(1-Ev*ev))*seb;
                 xc:=abp*(aa*xx+ff*yy) ;
                 yc:=abp*(bb*xx+gg*yy);
                 ro:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);
                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 teta:=tetag*180/pi;}
                 ti[1]:=abp;
                 ti[3]:=th;
                 ti[4]:=prd;

                 ti[5]:=aa;
                 ti[6]:=bb;
                 ti[7]:=ff;
                 ti[8]:=gg;

                 ti[9]:=omm;
                 ti[10]:=omb;
                 ti[11]:=ai;
                 ti[12]:=th;

                {writeln(otl2,ti[1]:25:10,' ',ti[2]:15:10,' ',ti[3]:15:10,' ',ti[4]:15:10,' ',
                ti[5]:15:10,' ',ti[6]:15:10,' ',ti[7]:15:10,' ',ti[8]:15:10);
                flush(otl2);}
                yc:=0;


end;

{procedure cirtpvd(ro0,teta0,rot,tetat,krvz,vr,massa,t,pia:extended; var ti:das);}
procedure cirtpvd(ro0,teta0,rot,tetat,krvz:extended;var vr,massa:extended;t,pia:extended; var ti:das);
var xs,ys,vx,vy,dxw,dyw,psi,mu,xa,ya,rt,za,xt,yt,zt,omm,omb,ai,e0,
comm,somm,comb,somb,ci,si,eb,sep,ceb,xx,yy,xc,yc,seb,tetag,rog,ev0,ro:extended;
 abp,ev,th,prd,aa,bb,ff,gg,omb2,zaf,cc,dd,a2,b2,c2,d2,aaa1,aaa2:extended;
 qq,m1,m2,vr1,vr2,z1,z2,vef,vef2,kp:extended;
begin
                kp:=0;
                xs:=ro0*sin(teta0*pi/180);
                ys:=ro0*cos(teta0*pi/180);
                vx:=rot*sin(teta0*pi/180)+ro0*tetat*cos(teta0*pi/180)*pi/180;
                vy:=rot*cos(teta0*pi/180)-ro0*tetat*sin(teta0*pi/180)*pi/180;
                teta0:=teta0*pi/180;

                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;

                mu:=sqrt(vx*vx+vy*vy);
                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                mu:=sqrt(vx*vx+vy*vy);
                {writeln(dt0,'mu=',mu:10:5,' psi=',psi*180/pi:10:5,' deg');}
                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                {kp:=(4*pi*pi)*massa;}
                qq:=ro0*krvz*abs(sin(psi-teta0))/(mu*mu);
                aa:=-(xs*vx+ys*vy)*(xs*vx+ys*vy);
                bb:=vx*vx+vy*vy;
                dd:=1/qq;
                cc:=-xs*xs-ys*ys;
                a2:=-dd;
                b2:=bb-dd*cc;
                c2:=bb*cc-aa;
                d2:=b2*b2-4*a2*c2;
                aaa1:=(-b2+sqrt(d2))/(2*a2);
                aaa2:=(-b2-sqrt(d2))/(2*a2);
                {kp=k*k}
                m1:=0;
                m2:=0;
                vr1:=0;
                vr2:=0;
                if aaa1>0 then
                begin
                     aaa1:=sqrt(aaa1);
                     m1:=aaa1*aaa1*aaa1/(qq*4*pi*pi);
                     z1:=aaa1*aaa1-xs*xs-ys*ys;
                     if z1>0 then
                     begin
                          z1:=sqrt(z1);
                          vr1:=(-xs*vx-ys*vy)/z1;
                     end else z1:=0;

                end;
                if aaa2>0 then
                begin
                     aaa2:=sqrt(aaa2);
                     m2:=aaa2*aaa2*aaa2/(qq*4*pi*pi);
                     kp:=aaa2*aaa2*aaa2/qq;
                     z2:=aaa2*aaa2-xs*xs-ys*ys;
                     if z2>0 then
                     begin
                          z2:=sqrt(z2);
                          vr2:=(-xs*vx-ys*vy)/z2;
                          vef:=vx*vx+vy*vy+vr2*vr2;
                          vef2:=kp/aaa2;
                          vef:=xs*xs+ys*ys+z2*z2;
                          aaa1:=0;
                     end;
                end;
                ti[1]:=m1;
                ti[2]:=vr1;
                ti[3]:=m2;
                ti[4]:=vr2;
                vr2:=vr2/((3600.0*24.0*365.2422)/149600000.0);
                {vr2:=vr2*1.00001;}
                bcir:=true;
                teta0:=teta0*180/pi;
                vr:=vr2;
                massa:=m2;

                {vr:=vr1;
                massa:=m1;}
                if krvz<0 then
                begin
                     m1:=0;
                end;
                if oprza then vr2:=-vr2;
                tpvd(ro0,teta0,rot,tetat,krvz,vr2,m2,0,1,ti);
                bcir:=false;
                m1:=0;

end;

(*procedure vydres;
var xs,ys,vx,vy,dxw,dyw,psi,mu,kp,xa,ya,rt,za,xt,yt,zt,abp,ev,omm,omb,ai,th,e0,prd,
comm,somm,comb,somb,ci,si,aa,bb,ff,gg,eb,sep,ceb,xx,yy,xc,yc,seb,tetag,rog,rok,xo,yo,mieph
,pex:extended;
dt0:text;
i,j,jp:integer;
begin
     assign(dt0,'okonchat.res');
     rewrite(dt0);
     {assign(rs2,'nev.res');
     rewrite(rs2);}
                j:=0;
                writeln(dt0,' idshodnye dannye:');
                writeln(dt0,'  rho0=',nro0:10:5,'+/-',sro0:10:5);
                writeln(dt0,'theta0=',nteta0:10:5,'+/-',steta0:10:5);
                writeln(dt0,'  rhot=',nrot:10:5,'+/-',srot:10:5);
                writeln(dt0,'thetat=',ntetat:10:5,'+/-',stetat:10:6);
                writeln(dt0);
                if bro0 then
                begin
                     j:=j+1;
                     writeln(dt0,'      rho0=',ro0:10:5,'+/-',sgm[j]:8:5,' var');
                end else
                begin
                     writeln(dt0,'      rho0=',ro0:10:5,'   not var');
                end;

                if bteta0 then
                begin
                     j:=j+1;
                     writeln(dt0,'    theta0=',teta0:8:3,'+/-',sgm[j]:6:3,' var');
                end else
                begin
                     writeln(dt0,'    theta0=',teta0:8:3,'   not var');
                end;

                if brot then
                begin
                     j:=j+1;
                     writeln(dt0,'      rhot=',rot:10:5,'+/-',sgm[j]:8:5,' var');
                end else
                begin
                     writeln(dt0,'      rhot=',rot:10:5,'   not var');
                end;

                if btetat then
                begin
                     j:=j+1;
                     writeln(dt0,'    thetat=',tetat:10:5,'+/-',sgm[j]:8:5,' var');
                end else
                begin
                     writeln(dt0,'    thetat=',tetat:10:5,'   not var');
                end;
                if br then
                begin
                     j:=j+1;
                     writeln(dt0,'  krivizna=',krvz:8:3,'+/-',sgm[j]:6:3,' var');
                end else
                begin
                     writeln(dt0,'  krivizna=',krvz:8:3,'   not var');

                end;
                if bvr then
                begin
                     j:=j+1;
                     writeln(dt0,'        vr=',vr:8:3,'+/-',sgm[j]:6:3,' var');
                end else
                begin
                     writeln(dt0,'        vr=',vr:8:3,'   not var');
                end;

                if bmassa then
                begin
                     j:=j+1;
                     writeln(dt0,'     massa=',massa:8:3,'+/-',sgm[j]:6:3,' var');
                end else
                begin
                     writeln(dt0,'     massa=',massa:8:3,'   not var');
                end;
                if bpia then
                begin
                     j:=j+1;
                     writeln(dt0,'       pia=',pia:8:3,'+/-',sgm[j]:6:3,' var');
                end else
                begin
                     writeln(dt0,'       pia=',pia:8:3,'   not var');
                end;

                {writeln(dt0,' krivizna=',krvz:10:5);
                writeln(dt0,'teta0=',(teta0):8:3);
                writeln(dt0,'   t0=',t0:8:3);
                writeln(dt0,'tetat=',tetat:20:15);
                writeln(dt0,'  rot=',rot:20:15);
                writeln(dt0,'   pi=',pia:8:3);
                writeln(dt0,'   vr=',vr:8:3);
                writeln(dt0,'massa=',massa:10:5);
                writeln(dt0);}

                xs:=ro0*sin(teta0*pi/180);
                ys:=ro0*cos(teta0*pi/180);
                vx:=rot*sin(teta0*pi/180)+ro0*tetat*cos(teta0*pi/180)*pi/180;
                vy:=rot*cos(teta0*pi/180)-ro0*tetat*sin(teta0*pi/180)*pi/180;
                teta0:=teta0*pi/180;

                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;

                mu:=sqrt(vx*vx+vy*vy);
                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                mu:=sqrt(vx*vx+vy*vy);
                writeln(dt0,'mu=',mu:10:5,' psi=',psi*180/pi:10:5,' deg');
                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                kp:=(4*pi*pi)*massa;
                xa:=xs/pia;
                ya:=ys/pia;
                rt:=por((kp*ro0*krvz*abs(sin(psi-teta0)))/(mu*mu),1/3);
                rok:=ro0/pia;
                writeln(dt0,'r=',rt:8:2,' a.e.  ro=',rok:8:2,' a.e.');
                writeln(dt0,' massa=',massa:4:1,' pi=',pia:8:4,'"');
                za:=sqrt(abs(rt*rt-rok*rok));
                if oprza then za:=-za;
                writeln(dt0,' za=',za:8:2,' a.e.');
                xt:=mu*sin(psi)/pia;
                yt:=mu*cos(psi)/pia;
                zt:=vr*(3600.0*24.0*365.2422)/149600000.0;
                kp:=sqrt(kp);
                elor(xa,ya,za,xt,yt,zt,kp,abp,ev,omm,omb,ai,th,e0);


                pvdv(xa,ya,za,xt,yt,zt,kp,abp,ev,ai,omb,omm,e0,th);

                writeln(dt0);
                if za<0 then writeln(dt0,' Z<0')
                        else writeln(dt0,' Z>0')  ;



                writeln(dt0,'a=',abp:10:5,'a.e.');
                writeln(dt0,'e=',ev:10:5);
                writeln(dt0,'omm=',omm*180/pi:10:5);
                writeln(dt0,'omb=',omb*180/pi:10:5);
                writeln(dt0,'i=',ai*180/pi:10:5);
                writeln(dt0,'T0=',(th+t0):10:5);
                writeln(dt0,'E0=',e0*180/pi:10:5);
                writeln(dt0,'a=',(abp*pia):10:5,'"');



                prd:=kp/(sqrt(abp*abp*abp)*2*pi);
                prd:=1/prd;
                abp:=abp*pia;
                writeln(dt0,'Period=',prd:10:5);
                writeln(dt0);
                writeln(dt0,' Ошибка ед веса по   ро=',oevr:10:5);
                writeln(dt0,' Ошибка ед веса по тета=',oevt:10:5);



            comm:=cos(omm);
            somm:=sin(omm);
            comb:=sin(omb);
            somb:=cos(omb);
            ci:=cos(ai);
            si:=sin(ai);
            aa:= comm*comb-somm*somb*ci;
            bb:= comm*somb+somm*comb*ci;
            ff:=-somm*comb-comm*somb*ci;
            gg:=-somm*somb+comm*comb*ci;
            for i:=1 to nt do
            begin
                 eb:=eqkepler((et[i]-th)*2*pi/prd,ev);
                 if (i=1) and fg then
                 begin
                      writeln('ev=',ev:10:5);
                      writeln('prd=',prd:10:5);
                      writeln('comm=',comm:10:5);
                      writeln('comb=',comb:10:5);
                      writeln('ci=',ci:10:5);
                      writeln('abp=',abp:10:5);
                      {halt(0);}
                 end;

                 seb:=sin(eb);
                 ceb:=cos(eb);
                 xx:=ceb-ev;
                 yy:=sqrt(1-Ev*ev)*seb;
                 xc:=abp*(aa*xx+ff*yy) ;
                 yc:=abp*(bb*xx+gg*yy);
                 rog:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);

                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 tetag:=tetag*180/pi;
                 xo:=ro[i]*sin(teta[i]*pi/180);
                 yo:=ro[i]*cos(teta[i]*pi/180);

                 writeln(dt0,i:4,' ',(et[i]+t0):8:3,' ',ro[i]:8:3,' ',teta[i]:8:3,' ',(ro[i]-rog):8:3,' ',(teta[i]-tetag):8:3);
                 {writeln(rs2,(et[i]+t0):8:3,' ',ro[i]:8:3,' ',teta[i]:8:3,' ',(ro[i]-rog):8:3,' ',
                 (teta[i]-tetag):8:3,' ',(xo-xc):8:3,' ',(yo-yc):8:3);}

            end;
            writeln(dt0);
            for i:=1 to niskl do
            begin
                 mieph:=1e39;

                 for j:=i to niskl do
                 begin
                      if ephi^[j]<mieph then
                      begin
                           mieph:=ephi^[j];
                           jp:=j;
                      end;
                 end;

                 pex:=ephi^[i];
                 ephi^[i]:=ephi^[jp];
                 ephi^[jp]:=pex;

                 pex:=roi^[i];
                 roi^[i]:=roi^[jp];
                 roi^[jp]:=pex;

                 pex:=tetai^[i];
                 tetai^[i]:=tetai^[jp];
                 tetai^[jp]:=pex;
            end;


            for i:=1 to niskl do
            begin
                 eb:=eqkepler((ephi^[i]-th)*2*pi/prd,ev);
                 seb:=sin(eb);
                 ceb:=cos(eb);
                 xx:=ceb-ev;
                 yy:=sqrt(1-Ev*ev)*seb;
                 xc:=abp*(aa*xx+ff*yy) ;
                 yc:=abp*(bb*xx+gg*yy);
                 rog:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);

                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 tetag:=tetag*180/pi;
                 writeln(dt0,(ephi^[i]+t0):8:3,' ',roi^[i]:8:3,' ',tetai^[i]:8:3,' ',(roi^[i]-rog):8:3,' ',
                 (tetai^[i]-tetag):8:3);

            end;
            writeln(dt0);

            if omm<0 then omm:=omm+2*pi;
            write(rzs,abp:10:5,' ',ev:10:5,' ',ai*180/pi:10:5,' ',omm*180/pi:10:5,' ',omb*180/pi:10:5,
                ' ',th+t0:10:5,' ',prd:10:5,' ',vr:10:5);

            for i:=2000 to t2050 do
            begin
                 eb:=eqkepler((i-t0-th)*2*pi/prd,ev);
                 seb:=sin(eb);
                 ceb:=cos(eb);
                 xx:=ceb-ev;
                 yy:=sqrt(1-Ev*ev)*seb;
                 xc:=abp*(aa*xx+ff*yy) ;
                 yc:=abp*(bb*xx+gg*yy);

                 rog:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);

                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 tetag:=tetag*180/pi;

                 write(rzs,' ',rog:8:3,' ',tetag:8:3);
            end;
            writeln(rzs);


            flush(dt0);
            close(dt0);

end;*)

procedure maxnev;
var xs,ys,vx,vy,dxw,dyw,psi,mu,kp,xa,ya,rt,za,xt,yt,zt,abp,ev,omm,omb,ai,th,e0,prd,
comm,somm,comb,somb,ci,si,aa,bb,ff,gg,eb,sep,ceb,xx,yy,xc,yc,seb,tetag,rog,rok,xo,yo
:extended;
i:integer;
begin


                xs:=ro0*sin(teta0*pi/180);
                ys:=ro0*cos(teta0*pi/180);
                vx:=rot*sin(teta0*pi/180)+ro0*tetat*cos(teta0*pi/180)*pi/180;
                vy:=rot*cos(teta0*pi/180)-ro0*tetat*sin(teta0*pi/180)*pi/180;
                teta0:=teta0*pi/180;

                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;

                mu:=sqrt(vx*vx+vy*vy);
                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                mu:=sqrt(vx*vx+vy*vy);
                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                kp:=(4*pi*pi)*massa;
                xa:=xs/pia;
                ya:=ys/pia;
                rt:=por((kp*ro0*krvz*abs(sin(psi-teta0)))/(mu*mu),1/3);
                rok:=ro0/pia;
                za:=sqrt(abs(rt*rt-rok*rok));
                if oprza then za:=-za;
                xt:=mu*sin(psi)/pia;
                yt:=mu*cos(psi)/pia;
                zt:=vr*(3600.0*24.0*365.2422)/149600000.0;
                kp:=sqrt(kp);
                elor(xa,ya,za,xt,yt,zt,kp,abp,ev,omm,omb,ai,th,e0);


                pvdv(xa,ya,za,xt,yt,zt,kp,abp,ev,ai,omb,omm,e0,th);




                prd:=kp/(sqrt(abp*abp*abp)*2*pi);
                prd:=1/prd;
                abp:=abp*pia;


            comm:=cos(omm);
            somm:=sin(omm);
            comb:=sin(omb);
            somb:=cos(omb);
            ci:=cos(ai);
            si:=sin(ai);
            aa:= comm*comb-somm*somb*ci;
            bb:= comm*somb+somm*comb*ci;
            ff:=-somm*comb-comm*somb*ci;
            gg:=-somm*somb+comm*comb*ci;
            oevr:=0;
            oevt:=0;
            for i:=1 to nt do
            begin
                 eb:=eqkepler((et[i]-th)*2*pi/prd,ev);
                 seb:=sin(eb);
                 ceb:=cos(eb);
                 xx:=ceb-ev;
                 yy:=sqrt(1-Ev*ev)*seb;
                 xc:=abp*(aa*xx+ff*yy) ;
                 yc:=abp*(bb*xx+gg*yy);
                 rog:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);

                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 tetag:=tetag*180/pi;
                 xo:=ro[i]*sin(teta[i]*pi/180);
                 yo:=ro[i]*cos(teta[i]*pi/180);
                 oevr:=oevr+sqr(rog-ro[i]);
                 oevt:=oevt+sqr(tetag-teta[i]);
                 {roo^[i]:=rog;
                 tetao^[i]:=tetag;}
            end;
            oevr:=sqrt(oevr/nt);
            oevt:=sqrt(oevt/nt);
            writeln('oevro=',oevr:10:5);
            writeln('oevte=',oevt:10:5);

end;


{procedure nevotap;
var xs,ys,vx,vy,dxw,dyw,psi,mu,kp,xa,ya,rt,za,xt,yt,zt,abp,ev,omm,omb,ai,th,e0,prd,
comm,somm,comb,somb,ci,si,aa,bb,ff,gg,eb,sep,ceb,xx,yy,xc,yc,seb,tetag,rog,rok,xo,yo:extended;
i,iop,j:integer;
    ron,tetan,yearn,nevron,nevten:array[1..700] of extended;
    apertn:array[1..700] of byte;
    chinap:array[0..255] of byte;
    sron,sten:array[0..255]  of extended;
    sro,ste:extended;
    dt,rs:text;

begin
                assign(dt,'mudr.res');
                reset(dt);
                iop:=0;
                while not(eof(dt)) do
                begin
                     iop:=iop+1;
                     readln(dt,yearn[iop],tetan[iop],ron[iop],apertn[iop]);
                end;
                close(dt);

                xs:=ro0*sin(teta0*pi/180);
                ys:=ro0*cos(teta0*pi/180);
                vx:=rot*sin(teta0*pi/180)+ro0*tetat*cos(teta0*pi/180)*pi/180;
                vy:=rot*cos(teta0*pi/180)-ro0*tetat*sin(teta0*pi/180)*pi/180;
                teta0:=teta0*pi/180;

                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;

                mu:=sqrt(vx*vx+vy*vy);
                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                mu:=sqrt(vx*vx+vy*vy);
                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                kp:=(4*pi*pi)*massa;
                xa:=xs/pia;
                ya:=ys/pia;
                rt:=por((kp*ro0*krvz*abs(sin(psi-teta0)))/(mu*mu),1/3);
                rok:=ro0/pia;
                za:=sqrt(abs(rt*rt-rok*rok));
                if oprza then za:=-za;
                xt:=mu*sin(psi)/pia;
                yt:=mu*cos(psi)/pia;
                zt:=vr*(3600.0*24.0*365.2422)/149600000.0;
                kp:=sqrt(kp);
                elor(xa,ya,za,xt,yt,zt,kp,abp,ev,omm,omb,ai,th,e0);
                pvdv(xa,ya,za,xt,yt,zt,kp,abp,ev,ai,omb,omm,e0,th);




                prd:=kp/(sqrt(abp*abp*abp)*2*pi);
                prd:=1/prd;
                abp:=abp*pia;


            comm:=cos(omm);
            somm:=sin(omm);
            comb:=sin(omb);
            somb:=cos(omb);
            ci:=cos(ai);
            si:=sin(ai);
            aa:= comm*comb-somm*somb*ci;
            bb:= comm*somb+somm*comb*ci;
            ff:=-somm*comb-comm*somb*ci;
            gg:=-somm*somb+comm*comb*ci;
            assign(rs,'kol.res');
            rewrite(rs);
            for i:=1 to iop do
            begin
                 eb:=eqkepler((yearn[i]-th-t0)*2*pi/prd,ev);
                 seb:=sin(eb);
                 ceb:=cos(eb);
                 xx:=ceb-ev;
                 yy:=sqrt(1-Ev*ev)*seb;
                 xc:=abp*(aa*xx+ff*yy) ;
                 yc:=abp*(bb*xx+gg*yy);
                 rog:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);

                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 tetag:=tetag*180/pi;

                 writeln(rs,yearn[i]:10:5,tetag:10:5,rog:10:5,tetan[i]:10:5,ron[i]:10:5);
                 nevron[i]:=rog-ron[i];
                 nevten[i]:=tetag-tetan[i];

            end;
            close(rs);
            for i:=0 to 255 do
            begin
                 chinap[i]:=0;
                 sron[i]:=0;
                 sten[i]:=0;
                 for j:=1 to iop do
                 begin
                      if i=apertn[j] then
                      begin
                           chinap[i]:=chinap[i]+1;
                           sron[i]:=sron[i]+sqr(nevron[j]);
                           sten[i]:=sten[i]+sqr(nevten[j]);
                      end;
                 end;
            end;
            assign(rs,'nap.res');
            rewrite(rs);
            for i:=0 to 255 do
            begin
                  if chinap[i]<>0 then
                  begin
                        sro:=sqrt(sron[i]/chinap[i]);
                        ste:=sqrt(sten[i]/chinap[i]);
                        writeln(rs,i:3,' ',sro:10:5,' ',ste:10:5,' ',chinap[i]);
                  end;
            end;
            close(rs);

end;}

procedure vydres2;
var xs,ys,vx,vy,dxw,dyw,psi,mu,kp,xa,ya,rt,za,xt,yt,zt,abp,ev,omm,omb,ai,th,e0,prd,
comm,somm,comb,somb,ci,si,aa,bb,ff,gg,eb,sep,ceb,xx,yy,xc,yc,seb,tetag,rog,rok,xo,yo
,tr0,sum,az,bz,hz,ad,bd,hd,etd:extended;
dt0,rs2,rs3:text;
i:integer;
begin
     assign(dt0,'nab.res');
     assign(rs2,'traekt.res');
     assign(rs3,'nev.res');
     rewrite(rs3);

     if brg then
     begin
                rewrite(dt0);
                rewrite(rs2);

                writeln(dt0,'   pi=',pia:8:3);
                {writeln(dt0,'massa=',massa:10:5);}
                writeln(dt0,'   t0=',t0:8:3);

                writeln(dt0);

                brg:=false;
                for i:=1  to nt do
                begin
                     xo:=ro[i]*sin(teta[i]*pi/180);
                     yo:=ro[i]*cos(teta[i]*pi/180);
                     writeln(rs2,et[i]+t0:8:3,' ',xo:8:3,' ',yo:8:3,' ',ro[i]:8:3,' ',teta[i]:8:3);
                end;
write(dt0,'massa  vr    kriv  theta0   ro0    mu      psi    rast    za      a    e       omm     omb      i      T0');
writeln(dt0,'    E0    period  oevr  oevt');

     end else
     begin
          append(dt0);
          append(rs2);
     end;



                {writeln(dt0,' krivizna=',krvz:10:5);
                writeln(dt0,'teta0=',(teta0):8:3);
                writeln(dt0,'  ro0=',ro0:8:3);
                writeln(dt0,'tetat=',tetat:20:15);
                writeln(dt0,'  rot=',rot:20:15);
                writeln(dt0);}
{8.619 317.562 6.590 0.03662  61.425 117.52  66.74 104.87 0.3019 217.613 228.949 145.385  1202 231.245 1073.9 0.1969}
                write(dt0,massa:4:1,' ',vr:6:3,' ',krvz:6:3,' ',teta0:7:3,' ',ro0:5:3);
                xs:=ro0*sin(teta0*pi/180);
                ys:=ro0*cos(teta0*pi/180);
                vx:=rot*sin(teta0*pi/180)+ro0*tetat*cos(teta0*pi/180)*pi/180;
                vy:=rot*cos(teta0*pi/180)-ro0*tetat*sin(teta0*pi/180)*pi/180;
                teta0:=teta0*pi/180;

                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;

                mu:=sqrt(vx*vx+vy*vy);
                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                mu:=sqrt(vx*vx+vy*vy);

                {writeln(dt0,'mu=',mu:10:5,' psi=',psi*180/pi:10:5,' deg');}
                write(dt0,' ',mu:7:5,' ',psi*180/pi:7:3);
                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                kp:=(4*pi*pi)*massa;
                xa:=xs/pia;
                ya:=ys/pia;
                rt:=por((kp*ro0*krvz*abs(sin(psi-teta0)))/(mu*mu),1/3);
                rok:=ro0/pia;

                {writeln(dt0,'r=',rt:8:2,' a.e.  ro=',rok:8:2,' a.e.');
                writeln(dt0,' massa=',massa:4:1,' pi=',pia:8:4,'"');}
                write(dt0,' ',rt:6:2);
                za:=sqrt(abs(rt*rt-rok*rok));
                if oprza then za:=-za;
                {writeln(dt0,' za=',za:8:2,' a.e.');}

                write(dt0,' ',za:6:2);

                xt:=mu*sin(psi)/pia;
                yt:=mu*cos(psi)/pia;
                zt:=vr*(3600.0*24.0*365.2422)/149600000.0;
                kp:=sqrt(kp);
                elor(xa,ya,za,xt,yt,zt,kp,abp,ev,omm,omb,ai,th,e0);
                pvdv(xa,ya,za,xt,yt,zt,kp,abp,ev,ai,omb,omm,e0,th);

                {writeln(dt0);
                if za<0 then writeln(dt0,' Z<0')
                        else writeln(dt0,' Z>0')  ;



                writeln(dt0,'a=',abp:10:5,'a.e.');
                writeln(dt0,'e=',ev:10:5);
                writeln(dt0,'omm=',omm*180/pi:10:5);
                writeln(dt0,'omb=',omb*180/pi:10:5);
                writeln(dt0,'i=',ai*180/pi:10:5);
                writeln(dt0,'T0=',(th+t0):10:5);
                writeln(dt0,'E0=',e0*180/pi:10:5);
                writeln(dt0,'a=',(abp*pia):10:5,'"');}
                write(dt0,' ',abp:6:2,' ',ev:6:4,' ',omm*180/pi:7:3,' ',omb*180/pi:7:3,' ',ai*180/pi:7:3,' '
                ,(th+t0):5:0,' ',e0*180/pi:7:3);



                prd:=kp/(sqrt(abp*abp*abp)*2*pi);
                prd:=1/prd;
                abp:=abp*pia;
                {writeln(dt0,'Period=',prd:10:5);}
                write(dt0,' ',prd:6:1,' ',oevr:8:6,' ',oevt:8:6);

            comm:=cos(omm);
            somm:=sin(omm);
            comb:=sin(omb);
            somb:=cos(omb);
            ci:=cos(ai);
            si:=sin(ai);
            aa:= comm*comb-somm*somb*ci;
            bb:= comm*somb+somm*comb*ci;
            ff:=-somm*comb-comm*somb*ci;
            gg:=-somm*somb+comm*comb*ci;
            az:=et[1]+t0;
            bz:=et[nt]+t0;
            hz:=(bz-az)/1000;
            for i:=0 to 1000 do
            begin
                 tr0:=az+i*hz;
                 {writeln(ev);}
                 eb:=eqkepler((tr0-th-t0)*2*pi/prd,ev);
                 seb:=sin(eb);
                 ceb:=cos(eb);
                 xx:=ceb-ev;
                 yy:=sqrt(1-Ev*ev)*seb;
                 xc:=abp*(aa*xx+ff*yy) ;
                 yc:=abp*(bb*xx+gg*yy);
                 rog:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);

                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 tetag:=tetag*180/pi;
                 {xo:=ro[i]*sin(teta[i]*pi/180);
                 yo:=ro[i]*cos(teta[i]*pi/180);}

                 writeln(rs2,tr0:8:3,' ',xc:20:10,' ',yc:20:10,' ',rog:20:10,' ',tetag:20:10);
                 {sum:=sum+sqr(xo-xc)+sqr(yo-yc);}
            end;
            sum:=0;
            for i:=1 to nt do
            begin
                 eb:=eqkepler((et[i]-th)*2*pi/prd,ev);
                 seb:=sin(eb);
                 ceb:=cos(eb);
                 xx:=ceb-ev;
                 yy:=sqrt(1-Ev*ev)*seb;
                 xc:=abp*(aa*xx+ff*yy) ;
                 yc:=abp*(bb*xx+gg*yy);
                 rog:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);

                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 tetag:=tetag*180/pi;
                 xo:=ro[i]*sin(teta[i]*pi/180);
                 yo:=ro[i]*cos(teta[i]*pi/180);
                 sum:=sum+sqr(xo-xc)+sqr(yo-yc);
                 {writeln(rs3,et[i]+t0:8:3,' ',ro[i]:8:3,' ', teta[i]:9:3,' ',xo:8:3,' ',yo:8:3);}
            end;
            sum:=sqrt(sum/(i*2));
            writeln(dt0,sum:10:5);
            ad:=th;
            bd:=th+prd;
            hd:=(bd-ad)/1000;

            for i:=0 to 1000 do
            begin
                 etd:=ad+hd*i;
                 eb:=eqkepler((etd-th)*2*pi/prd,ev);
                 seb:=sin(eb);
                 ceb:=cos(eb);
                 xx:=ceb-ev;
                 yy:=sqrt(1-Ev*ev)*seb;
                 xc:=abp*(aa*xx+ff*yy) ;
                 yc:=abp*(bb*xx+gg*yy);
                 rog:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);

                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 tetag:=tetag*180/pi;
                 writeln(rs3,etd+t0:8:3,' ',rog:8:3,' ', tetag:9:3,' ',xc:8:3,' ',yc:8:3);
            end;
            close(rs3);
            {teta0:=teta0*180/pi;
            for i:=1 to nt do
            begin
                 tpvd(ro0,teta0,rot,tetat,krvz,vr,massa,et[i],pia, rog,tetag,t0 );

                 writeln(dt0,(et[i]+t0):8:3,' ',ro[i]:8:3,' ',teta[i]:8:3,' ',(ro[i]-rog):8:3,' ',(teta[i]-tetag):8:3);

            end;}
            close(rs2);
            flush(dt0);
            close(dt0);

end;


procedure gen0;
var xs,ys,vx,vy,dxw,dyw,psi,mu,kp,xa,ya,rt,za,xt,yt,zt,abp,ev,omm,omb,ai,th,e0,prd,
comm,somm,comb,somb,ci,si,aa,bb,ff,gg,eb,sep,ceb,xx,yy,xc,yc,seb,tetag,rog,rok,xo,yo
,tr0,sum,az,bz,hz:extended;
dt0,rs2:text;
i:integer;
begin

                xs:=ro0*sin(teta0*pi/180);
                ys:=ro0*cos(teta0*pi/180);
                vx:=rot*sin(teta0*pi/180)+ro0*tetat*cos(teta0*pi/180)*pi/180;
                vy:=rot*cos(teta0*pi/180)-ro0*tetat*sin(teta0*pi/180)*pi/180;
                teta0:=teta0*pi/180;

                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;

                mu:=sqrt(vx*vx+vy*vy);
                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                mu:=sqrt(vx*vx+vy*vy);

                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                kp:=(4*pi*pi)*massa;
                xa:=xs/pia;
                ya:=ys/pia;
                rt:=por((kp*ro0*krvz*abs(sin(psi-teta0)))/(mu*mu),1/3);
                rok:=ro0/pia;

                za:=sqrt(abs(rt*rt-rok*rok));
                if oprza then za:=-za;


                xt:=mu*sin(psi)/pia;
                yt:=mu*cos(psi)/pia;
                zt:=vr*(3600.0*24.0*365.2422)/149600000.0;
                kp:=sqrt(kp);
                elor(xa,ya,za,xt,yt,zt,kp,abp,ev,omm,omb,ai,th,e0);
                pvdv(xa,ya,za,xt,yt,zt,kp,abp,ev,ai,omb,omm,e0,th);
                prd:=kp/(sqrt(abp*abp*abp)*2*pi);
                prd:=1/prd;
                abp:=abp*pia;
                if omm<0 then omm:=omm+2*pi;

                abp0:=abp;
                ev0:=ev;
                ai0:=ai;
                omm0:=omm;

                omb0:=omb;
                th0:=th;
                prd0:=prd;



            comm:=cos(omm);
            somm:=sin(omm);
            comb:=sin(omb);
            somb:=cos(omb);
            ci:=cos(ai);
            si:=sin(ai);
            aa:= comm*comb-somm*somb*ci;
            bb:= comm*somb+somm*comb*ci;
            ff:=-somm*comb-comm*somb*ci;
            gg:=-somm*somb+comm*comb*ci;
            for i:=1 to nt do
            begin
                 eb:=eqkepler((et[i]-th)*2*pi/prd,ev);
                 seb:=sin(eb);
                 ceb:=cos(eb);
                 xx:=ceb-ev;
                 yy:=sqrt(1-Ev*ev)*seb;
                 xc:=abp*(aa*xx+ff*yy) ;
                 yc:=abp*(bb*xx+gg*yy);
                 xro0^[i]:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);

                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 xteta0^[i]:=tetag*180/pi;
                 {writeln(otl, et[i]:10:3,' ',xro0^[i]:7:3,' ',xteta0^[i]:8:3);}

            end;

            for i:=2000 to t2050 do
            begin
                 eb:=eqkepler((i-t0-th)*2*pi/prd,ev);
                 seb:=sin(eb);
                 ceb:=cos(eb);
                 xx:=ceb-ev;
                 yy:=sqrt(1-Ev*ev)*seb;
                 xc:=abp*(aa*xx+ff*yy) ;
                 yc:=abp*(bb*xx+gg*yy);
                 roe^[i]:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);

                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 tetae^[i]:=tetag*180/pi;
                 {writeln(otl, et[i]:10:3,' ',xro0^[i]:7:3,' ',xteta0^[i]:8:3);}

            end;



end;

procedure anres;
var xs,ys,vx,vy,dxw,dyw,psi,mu,kp,xa,ya,rt,za,xt,yt,zt,abp,omm,omb,ai,th,e0,prd,
comm,somm,comb,somb,ci,si,aa,bb,ff,gg,eb,sep,ceb,xx,yy,xc,yc,seb,tetag,rog,rok,xo,yo
,tr0,sum,az,bz,hz:extended;
dt0,rs2:text;
i:integer;
begin
                if (ev>1) or (krvz<0.1) then
                begin
                     iv:=iv+1;
                     exit;
                end;
                xs:=ro0*sin(teta0*pi/180);
                ys:=ro0*cos(teta0*pi/180);
                vx:=rot*sin(teta0*pi/180)+ro0*tetat*cos(teta0*pi/180)*pi/180;
                vy:=rot*cos(teta0*pi/180)-ro0*tetat*sin(teta0*pi/180)*pi/180;
                teta0:=teta0*pi/180;

                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;

                mu:=sqrt(vx*vx+vy*vy);
                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                mu:=sqrt(vx*vx+vy*vy);

                dxw:=vx;
                dyw:=vy;
                psi:=arctan(dxw/dyw);
                if dyw<0 then psi:=pi+psi;
                if psi<0 then psi:=psi+2*pi;
                kp:=(4*pi*pi)*massa;
                xa:=xs/pia;
                ya:=ys/pia;
                rt:=por((kp*ro0*krvz*abs(sin(psi-teta0)))/(mu*mu),1/3);
                rok:=ro0/pia;

                za:=sqrt(abs(rt*rt-rok*rok));
                if oprza then za:=-za;


                xt:=mu*sin(psi)/pia;
                yt:=mu*cos(psi)/pia;
                zt:=vr*(3600.0*24.0*365.2422)/149600000.0;
                kp:=sqrt(kp);
                elor(xa,ya,za,xt,yt,zt,kp,abp,ev,omm,omb,ai,th,e0);
                pvdv(xa,ya,za,xt,yt,zt,kp,abp,ev,ai,omb,omm,e0,th);
                prd:=kp/(sqrt(abp*abp*abp)*2*pi);
                prd:=1/prd;
                abp:=abp*pia;

                if omm<0 then omm:=omm+2*pi;

                sabp0:=sabp0+sqr(abp0-abp);
                sev0:=sev0+sqr(ev0-ev);
                sai0:=sai0+sqr(ai0-ai);
                {if (omm0-omm)>pi then omm:=omm+pi;}
                somm0:=somm0+sqr(omm0-omm);
                somb0:=somb0+sqr(omb0-omb);
                sth0:=sth0+sqr(th0-th);
                sprd0:=sprd0+sqr(prd0-prd);


                writeln(otl,abp:10:5,' ',ev:10:5,' ',ai*180/pi:10:5,' ',omm*180/pi:10:5,' ',omb*180/pi:10:5,
                ' ',th:10:5,' ',prd:10:5,' ',vr:10:5);
                flush(otl);
                write(rzs,abp:10:5,' ',ev:10:5,' ',ai*180/pi:10:5,' ',omm*180/pi:10:5,' ',omb*180/pi:10:5,
                ' ',th+t0:10:5,' ',prd:10:5,' ',vr:10:5);



            comm:=cos(omm);
            somm:=sin(omm);
            comb:=sin(omb);
            somb:=cos(omb);
            ci:=cos(ai);
            si:=sin(ai);
            aa:= comm*comb-somm*somb*ci;
            bb:= comm*somb+somm*comb*ci;
            ff:=-somm*comb-comm*somb*ci;
            gg:=-somm*somb+comm*comb*ci;

            for i:=2000 to t2050 do
            begin
                 eb:=eqkepler((i-t0-th)*2*pi/prd,ev);
                 seb:=sin(eb);
                 ceb:=cos(eb);
                 xx:=ceb-ev;
                 yy:=sqrt(1-Ev*ev)*seb;
                 xc:=abp*(aa*xx+ff*yy) ;
                 yc:=abp*(bb*xx+gg*yy);

                 rog:=sqrt(xc*xc+yc*yc);

                 tetag:=arctan(xc/yc);

                 if yc<0 then tetag:=pi+tetag;
                 if tetag<0 then tetag:=tetag+2*pi;
                 tetag:=tetag*180/pi;

                 sroe^[i]:=sroe^[i]+sqr(roe^[i]-rog);
                 statae^[i]:=statae^[i]+sqr(tetae^[i]-tetag);
                 write(rzs,' ',rog:8:3,' ',tetag:8:3);
            end;
            writeln(rzs,' ',massa:8:4);


end;
procedure vydres3(t0,pia:extended);
var dt0:text;
begin
     close(rzs);
     writeln(otl,' procedure vydres3');
     flush(otl);
     sabp0:=sqrt(sabp0/(t100-iv));
     sev0:=sqrt(sev0/(t100-iv));
     sai0:=sqrt(sai0/(t100-iv));
     somm0:=sqrt(somm0/(t100-iv));
     somb0:=sqrt(somb0/(t100-iv));
     sth0:=sqrt(sth0/(t100-iv));
     sprd0:=sqrt(sprd0/(t100-iv));
     writeln(otl,'ku');
     flush(otl);

     assign(dt0,'okonchat.res');
     append(dt0);
     writeln(dt0);
     writeln(dt0,'a=   ',abp0:10:5,'"');
     writeln(dt0,'  +/-',sabp0:10:5,'"');
     writeln(dt0);


     writeln(dt0,'e=    ',ev0:10:5);
     writeln(dt0,'   +/-',sev0:10:5);
     writeln(dt0);

     writeln(dt0,'omm=   ',omm0*180/pi:10:5);
     writeln(dt0,'    +/-',somm0*180/pi:10:5);
     writeln(dt0);

     writeln(dt0,'omb=   ',omb0*180/pi:10:5);
     writeln(dt0,'    +/-',somb0*180/pi:10:5);
     writeln(dt0);


     writeln(dt0,'i=   ',ai0*180/pi:10:5);
     writeln(dt0,'  +/-',sai0*180/pi:10:5);
     writeln(dt0);


     writeln(dt0,'T0=   ',(th0+t0):10:5);
     writeln(dt0,'   +/-',(sth0):10:5);
     writeln(dt0);

     writeln(dt0,'prd=   ',prd0:10:5);
     writeln(dt0,'    +/-',sprd0:10:5);
     writeln(dt0);



     writeln(dt0);
     writeln(dt0);
     writeln(otl,'ku2');
     flush(otl);

     writeln(dt0,'Ефмериды:');
     writeln(dt0);
     writeln(dt0);
     for i:=2000 to t2050 do
     begin
           sroe^[i]:=sqrt(sroe^[i]/(t100-iv));
           statae^[i]:=sqrt(statae^[i]/(t100-iv));
           writeln(dt0,i:4,'.0 ',roe^[i]:7:3,' ',sroe^[i]:7:3,' ',tetae^[i]:8:3,' ',statae^[i]:7:3);
     end;
     close(dt0);
     writeln(otl,'ku3');
     flush(otl);


end;

procedure vgal1(var a,e,omm,omb,ai,v,akapa, x,y,z,xt,yt,zt:extended);

var u, sqp,r,
    rt,rvt,aa,b,c:extended;
begin
     { Procedure для вычисления относительного положения и скорости };
      {   двойной звезды по элементам орбиты. Используются procedures rot и rot2 };
      { Часть алгоритм взята из 'Введения в теор. астрономию' Субботина стр. 78 };
      {  x,y,z радиус вектор тела                 };
      {  xt,yt,zt вектор скорости                 };
      {  v истинная аномалия                      };
      {  akapa=sqrt(m*g) для относительного движения    };
      {  a большая полуось                        };
      {  e экцентриситет                          };
      {  omm аргумент перицентра                  };
      {  omb долгота восходящего узла             };
      {  ai наклон орбиты                         };
     { assign(lv,'izty.pas');
      reset(lv);
      readln(lv,a,e,omm,omb,ai,v,amt,g );     }
      r:=a*(1-e*e)/(1+e*cos(v)) ;
      sqp:=sqrt(a*(1-e*e));
      rt:=akapa*e*sin(v)/sqp;
      rvt:=akapa*sqp/r;
      u:=omm+v;
      aa:=(cos(u)*cos(omb)-sin(u)*sin(omb)*cos(ai));
      b:=(cos(u)*sin(omb)+sin(u)*cos(omb)*cos(ai));
      c:=sin(u)*sin(ai) ;
      x:=r*aa;
      y:=r*b;
      z:=r*c;
      xt:=rt*aa+rvt*(-sin(u)*cos(omb)-cos(u)*sin(omb)*cos(ai));
      yt:=rt*b+rvt*(-sin(u)*sin(omb)+cos(u)*cos(omb)*cos(ai));
      zt:=rt*c+rvt*sin(ai)*cos(u);
end;
procedure GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,t2:extended;var iro0,iteta0,rot,tetat,r,vr,massa:extended);
var i:integer;
v,akapa,x,y,z,xt,yt,zt,mm,eb,ai3,eps:extended;
abp2,ev2,omm2,omb2,ai2,th2:extended;
x2,y2,ebc,cebc,sebc,xxc,yyc,omb3,x3,y3,z3,z2,t3,vx,vy:extended;
dxw,dyw,ro2,teta2,ro3,teta3,psi,r3,mu:extended;

begin
     eps:=1e-8;

     massa:=abp*abp*abp/(prd*prd);

     mm:=(t2-th)*2*pi/prd;
     eb:=eqkepler(mm,ev);
     v:=sqrt(ABS((1+ev)/(1-ev)))*tg(eb/2);
     v:=arctan(v);
     v:=2*v;
     akapa:=((4*pi*pi)*massa);

     omb3:=pi/2-omb;
     vgal1(abp,ev,omm,omb3,ai,v,akapa,x,y,z,xt,yt,zt);

     dxw:=x;
     dyw:=y;
     iro0:=sqrt(sqr(dyw)+sqr((dxw)));
     iteta0:=arctan(dxw/dyw);
     if dyw<0 then iteta0:=pi+iteta0;
     if iteta0<0 then iteta0:=iteta0+2*pi;
     iteta0:=iteta0*180/pi;

     {rtti2(ti0,t2,x2,y2,ebc,cebc,sebc,xxc,yyc);}

     t3:=t2-eps;
     mm:=(t3-th)*2*pi/prd;
     eb:=eqkepler(mm,ev);
     v:=sqrt((1+ev)/(1-ev))*tg(eb/2);
     v:=arctan(v);
     v:=2*v;

     vgal1(abp,ev,omm,omb3,ai,v,akapa,x2,y2,z2,xt,yt,zt);
     dxw:=x2;
     dyw:=y2;
     ro2:=sqrt(sqr(dyw)+sqr((dxw)));
     teta2:=arctan(dxw/dyw);
     if dyw<0 then teta2:=pi+teta2;
     if teta2<0 then teta2:=teta2+2*pi;
     teta2:=teta2*180/pi;

     t3:=t2+eps;
     mm:=(t3-th)*2*pi/prd;
     eb:=eqkepler(mm,ev);
     v:=sqrt((1+ev)/(1-ev))*tg(eb/2);
     v:=arctan(v);
     v:=2*v;

     vgal1(abp,ev,omm,omb3,ai,v,akapa,x3,y3,z3,xt,yt,zt);
     dxw:=x3;
     dyw:=y3;
     ro3:=sqrt(sqr(dyw)+sqr((dxw)));
     teta3:=arctan(dxw/dyw);
     if dyw<0 then teta3:=pi+teta3;
     if teta3<0 then teta3:=teta3+2*pi;
     teta3:=teta3*180/pi;

     zt:=(z3-z2)/(2*eps);
     vr:=zt/((3600.0*24.0*365.2422)/149600000.0);
     rot:=(ro3-ro2)/(2*eps);
     tetat:=(teta3-teta2)/(2*eps);

     vx:=rot*sin(iteta0*pi/180)+iro0*tetat*cos(iteta0*pi/180)*pi/180;
     vy:=rot*cos(iteta0*pi/180)-iro0*tetat*sin(iteta0*pi/180)*pi/180;

     dxw:=vx;
     dyw:=vy;
     psi:=arctan(dxw/dyw);
     if dyw<0 then psi:=pi+psi;
     if psi<0 then psi:=psi+2*pi;
     mu:=(vx*vx+vy*vy);
     r3:=sqrt(x*x+y*y+z*z);
     r3:=r3/pia;
     r3:=r3*r3*r3;
     r:=(akapa*iro0*sin(psi-iteta0*pi/180));
     if r<>0 then
     begin
          r:=abs(r3*mu/r);
     end else
     begin
          r:=10000;
     end;

end;

procedure DrgEphPvd(iro0,iteta0,rot,tetat,r,vr,massa,t2:extended;
                    var diro0,diteta0,drot,dtetat,dr,dvr,dmassa:extended);
var ti0,ti:das;
i:integer;
abp,e,omm,omb,ai,v,akapa,x,y,z,xt,yt,zt,mm,eb,prd,ai3,eps:extended;
abp2,ev2,omm2,omb2,ai2,th2:extended;
x2,y2,ebc,cebc,sebc,xxc,yyc,omb3,x3,y3,z3,z2,t3,vx,vy:extended;
dxw,dyw,ro2,teta2,ro3,teta3,psi,r3,mu:extended;
begin
      eps:=1e-8;
      if mbcir then
      begin
           cirtpvd(iro0,iteta0,rot,tetat,r,vr,massa,0,pia,ti0);
      end else
      begin
           tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);
      end;
     abp:=ti0[1];
     e:= ti0[2];
     th:= ti0[3];
     prd:=ti0[4];
     omm:=ti0[9];
     omb:=ti0[10];
     ai:= ti0[11];
     ev:=e;
     {t2:=0;}



     mm:=(t2-th)*2*pi/prd;
     eb:=eqkepler(mm,ev);
     v:=sqrt((1+e)/(1-e))*tg(eb/2);
     v:=arctan(v);
     v:=2*v;
     akapa:=((4*pi*pi)*massa);

     omb3:=pi/2-omb;
     vgal1(abp,e,omm,omb3,ai,v,akapa,x,y,z,xt,yt,zt);

     dxw:=x;
     dyw:=y;
     diro0:=sqrt(sqr(dyw)+sqr((dxw)));
     diteta0:=arctan(dxw/dyw);
     if dyw<0 then diteta0:=pi+diteta0;
     if diteta0<0 then diteta0:=diteta0+2*pi;
     diteta0:=diteta0*180/pi;

     {rtti2(ti0,t2,x2,y2,ebc,cebc,sebc,xxc,yyc);}

     t3:=t2-eps;
     mm:=(t3-th)*2*pi/prd;
     eb:=eqkepler(mm,ev);
     v:=sqrt((1+e)/(1-e))*tg(eb/2);
     v:=arctan(v);
     v:=2*v;

     vgal1(abp,e,omm,omb3,ai,v,akapa,x2,y2,z2,xt,yt,zt);
     dxw:=x2;
     dyw:=y2;
     ro2:=sqrt(sqr(dyw)+sqr((dxw)));
     teta2:=arctan(dxw/dyw);
     if dyw<0 then teta2:=pi+teta2;
     if teta2<0 then teta2:=teta2+2*pi;
     teta2:=teta2*180/pi;

     t3:=t2+eps;
     mm:=(t3-th)*2*pi/prd;
     eb:=eqkepler(mm,ev);
     v:=sqrt((1+e)/(1-e))*tg(eb/2);
     v:=arctan(v);
     v:=2*v;

     vgal1(abp,e,omm,omb3,ai,v,akapa,x3,y3,z3,xt,yt,zt);
     dxw:=x3;
     dyw:=y3;
     ro3:=sqrt(sqr(dyw)+sqr((dxw)));
     teta3:=arctan(dxw/dyw);
     if dyw<0 then teta3:=pi+teta3;
     if teta3<0 then teta3:=teta3+2*pi;
     teta3:=teta3*180/pi;

     zt:=(z3-z2)/(2*eps);
     dvr:=zt/((3600.0*24.0*365.2422)/149600000.0);
     drot:=(ro3-ro2)/(2*eps);
     dtetat:=(teta3-teta2)/(2*eps);

     vx:=drot*sin(diteta0*pi/180)+diro0*dtetat*cos(diteta0*pi/180)*pi/180;
     vy:=drot*cos(diteta0*pi/180)-diro0*dtetat*sin(diteta0*pi/180)*pi/180;

     dxw:=vx;
     dyw:=vy;
     psi:=arctan(dxw/dyw);
     if dyw<0 then psi:=pi+psi;
     if psi<0 then psi:=psi+2*pi;
     mu:=(vx*vx+vy*vy);
     r3:=sqrt(x*x+y*y+z*z);
     r3:=r3/pia;
     r3:=r3*r3*r3;
     dr:=abs(r3*mu/(akapa*diro0*sin(psi-diteta0*pi/180)));
     {tpvd(diro0,diteta0,drot,dtetat,dr,dvr,massa,et^[1],pia,ti);
     i:=0;}
end;
procedure XyzFrmEl(a,v,e,ai,omb,omm:extended;var x,y,z:extended);
var u, sqp,r,
    aa,b,c:extended;
begin
      r:=a*(1-e*e)/(1+e*cos(v)) ;
      sqp:=sqrt(a*(1-e*e));
      u:=omm+v;
      aa:=(cos(u)*cos(omb)-sin(u)*sin(omb)*cos(ai));
      b:=(cos(u)*sin(omb)+sin(u)*cos(omb)*cos(ai));
      c:=sin(u)*sin(ai) ;
      x:=r*aa;
      y:=r*b;
      z:=r*c;

end;

function GetUgol(abp1,ev1,th1,prd1,ai1,omb1,omm1,
                  abp2,ev2,th2,prd2,ai2,omm2,omb2:extended):extended;
var ax,ay,az,bx,by,bz,v,x1,y1,z1,x2,y2,z2,sk,m1,m2,cosug:extended;

begin
      v:=0;
      XyzFrmEl(abp1,v,ev1,ai1,omb1,omm1,ax,ay,az);
      v:=pi/2;
      XyzFrmEl(abp1,v,ev1,ai1,omb1,omm1,bx,by,bz);
      x1:=ay*bz-az*by;
      y1:=az*bx-ax*bz;
      z1:=ax*by-ay*bx;
      v:=0;
      XyzFrmEl(abp2,v,ev2,ai2,omb2,omm2,ax,ay,az);
      v:=pi/2;
      XyzFrmEl(abp2,v,ev2,ai2,omb2,omm2,bx,by,bz);
      x2:=ay*bz-az*by;
      y2:=az*bx-ax*bz;
      z2:=ax*by-ay*bx;
      sk:=x1*x2+y1*y2+z1*z2;
      m1:=sqrt(x1*x1+y1*y1+z1*z1);
      m2:=sqrt(x2*x2+y2*y2+z2*z2);
      cosug:=sk/(m1*m2);
      GetUgol:=arccos(cosug);

end;
procedure ogr2(var raz:extended);
var i,mi:integer;
min,s:extended;
begin
     mi:=0;
     min:=1e39;
     for i:=-10 to 10 do
     begin
          s:=abs(raz+i*360);
          if s<min then
          begin
               min:=s;
               mi:=i;
          end;
     end;
     raz:=raz+mi*360;
end;

procedure  vydnev(tb:integer);
var tvr2,tmassa,dro,dte,x,y,xg,yg:extended;
ti0:das;
ivyh:integer;
begin
      if mbcir then
      begin
           mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,tvr2,tmassa,et^[1],1,ti0);
      end else
      begin
           mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);
      end;
      if not(nevyh) then ivyh:=1 else ivyh:=0;

      for id:=1 to nt do
      begin
           {if gbop^[id]<>0 then}
           begin
                mpdq.rtti(ti0,et^[id],rog,tetag);

                dro:=ro^[id]-rog;
                dte:=(teta^[id]-tetag);
                ogr2(dte);
                x:=ro^[id]*sin(teta^[id]*pi/180);
                y:=ro^[id]*cos(teta^[id]*pi/180);
                xg:=rog*sin(tetag*pi/180);
                yg:=rog*cos(tetag*pi/180);
                ogr2(DTE);

                write(arsv[tb],aref^[id],' ',ivyh:1,' ',swds,' ',comp10,' ',m1:7:3,' ',m2:7:3
                ,' ',et^[id]+t0:8:3,' ',gbop^[id]:1,' ',gbop3^[id]:1,' '
                ,ro^[id]:8:3,' ',teta^[id]:8:3,' ',rog:8:3,' ',tetag:8:3,' '
               ,dro:8:3,' ',dte:8:3,' ',x:8:3,' ',y:8:3,' ',gf0:15:10,' ',ti0[2]:15:10);
               if btvht then
               begin
                    write(arsv[tb],' 1');

               end else
               begin
                    write(arsv[tb],' 0');

               end;
               writeln(arsv[tb],' ',l);

               {clrscr;
               writeln(l);
               readln;}
           end;

      end;
      flush(arsv[tb]);


end;
function GetUgol2(abp1,ev1,th1,prd1,ai1,omb1,omm1,
                          abp2,ev2,th2,prd2,ai2,omm2,omb2:extended):extended;
var ug1,ug2:extended;
begin
     ug1:=GetUgol(abp1,ev1,th1,prd1,ai1,omb1,omm1,
                          abp2,ev2,th2,prd2,ai2,omm2,omb2);
     if ug1>pi/2 then ug1:=pi-ug1;

     ug2:=GetUgol(abp1,ev1,th1,prd1,ai1,omb1+pi,omm1+pi,
                          abp2,ev2,th2,prd2,ai2,omm2,omb2);
     if ug2>pi/2 then ug2:=pi-ug2;
     if ug2<ug1 then ug1:=ug2;
     GetUgol2:=ug1;
end;

begin
     bcir:=false;
     fg:=false;
     brg:=true;
     sabp0:=0;
     sev0:=0;
     sai0:=0;
     somm0:=0;
     somb0:=0;
     sth0:=0;
     sprd0:=0;
     new(sroe); new(statae);
     new(roe); new(tetae);
     {new(roi); new(tetai); new(ephi);}
     for jh:=2000 to t2050 do sroe^[jh]:=0;
     for jh:=2000 to t2050 do statae^[jh]:=0;
     iv:=0;
     {assign(rzs,'elorvr.res');
     rewrite(rzs);

     assign(otl2,'mpdq.otl');
     rewrite(otl2);}

end.
