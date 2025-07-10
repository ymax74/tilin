{$A+,B-,D+,E-,F+,G+,I+,L+,N+,O+,P-,Q+,R+,S+,T-,V+,X+,Y+}
{$M 65520,0,655360}
uses sbop,derty,
mnk,mpvd,mpdq,crt,dos,ubrk,derty2,getelm;
const ovds2='..\ovdstr2';
var svpst:string;



{procedure ogr1(var raz:extended);
begin
     if raz<-180 then raz:=raz+360;
     if raz> 180 then raz:=raz+360;
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

procedure erval(psd:string; code:integer);
begin
     if code<>0 then
     begin
          writeln(' Error in procedure val.');
          writeln('string=',psd);
          readln;
     end;
end;}




function otnc(rf:string):string;
var pst:string;
   k,i:integer;
   th,tc,tm:boolean;
begin
     k:=length(rf);
     i:=0;
     th:=true;
     pst[0]:=chr(k);
     tm:=false;
     while (i<=k) and th do
     begin
          i:=i+1;
          tc:= (rf[i]='-') or (rf[i]='1') or (rf[i]='2') or (rf[i]='3') or (rf[i]='4') or (rf[i]='5') or
          (rf[i]='6') or (rf[i]='7') or (rf[i]='8') or (rf[i]='9') or  (rf[i]='0') or (rf[i]='.');
          if tc then pst[i]:=rf[i]
                else pst[i]:=' ';
          if tc then tm:=true;
          if not(tc) and tm then th:=false;
     end;
     pst[0]:=chr(i-1);
     otnc:=pst;
end;



{function eqkepler(m,e:extended):extended;
const eps=1e-8;
var ep,eb,ep0:extended;
begin
      ep:=m;
      repeat
            eb:=m+e*sin(ep);
            ep0:=ep;
            ep:=eb;
      until  abs(ep0-eb)<eps;
      eqkepler:=eb;
end;}
function por(x,p:extended):extended;
begin
     por:=exp(p*ln(x));
end;

(*procedure iskl;
const t3=10;
var mxt,mxr:extended;
    imxr,imxt,i:integer;
begin
     mxt:=0;
     mxr:=0;

     for i:=1 to nt do
     begin
          if abs(roo^[i]-ro^[i])>mxr then
          begin
               mxr:=abs(roo^[i]-ro^[i]);
               imxr:=i;
          end;
          if abs(tetao^[i]-teta^[i])>mxt then
          begin
               mxt:=abs(tetao^[i]-teta^[i]);
               imxt:=i;
          end;
     end;
     mxr:=mxr/oevr;
     mxt:=mxt/oevt;
     if mxt>mxr then
     begin
          mxr:=mxt;
          imxr:=imxt;
     end;
     ew:=mxr>t3;

     if ew then cvb(imxr);


end;*)


function ctg(x:extended):extended;
begin
     ctg:=cos(x)/sin(x);
end;
{procedure erval(fed:string; ercode:integer);
begin
     if ercode>0 then
     begin
          writeln(' stroka=',fed,' code=',ercode);
          readln;
     end;
end;}

function cirnev(r2,om,om0,x0,y0:extended):extended;
var km:integer;
xc,yc,nevx,nevy,smr,sn,cs:extended;
begin
     km:=0;
     smr:=0;
           for i:=1 to nnt do
           begin
                if gbop^[i]=0 then
                begin
                     km:=km+1;
                     sn:=sin(om*et^[i]+om0);
                     cs:=cos(om*et^[i]+om0);
                     xc:=r2*sn+x0;
                     yc:=r2*cs+y0;
                     nevx:=x^[i]-xc;
                     nevy:=y^[i]-yc;
                     smr:=smr+sqr(nevx)+sqr(nevy);
                end;
           end;
           cirnev:=smr;

end;
function cirlam(l,r2,om,om0,x0,y0,dr2,dom,dom0,dx0,dy0:extended):extended;
begin
     r2 := r2-l*dr2;
     if r2>0 then
     begin
          om := om-l*dom;
          om0:=om0-l*dom0;
          x0 := x0-l*dx0;
          y0 := y0-l*dy0;
          cirlam:=cirnev(r2,om,om0,x0,y0);
     end else
     begin
           cirlam:=1e39;
     end;

end;

(*procedure flcir(var r2,om,om0,x0,y0:extended);
var l,xc,yc,nevx,nevy,smr,sn,cs:extended;
rs:text;
itr:integer;
vyh:boolean;
br2,bom,bom0,bx0,by0,sml,smn,raz:extended;
nsl:boolean;
label l1;


begin
     if r2>3000 then
     begin
           itrf:=3001;
           exit;
     end;
     l:=0.05;
     {om:=-mu/r;}
     {om0:=teta0*pi/180;}
     {X0:=-X0;
     y0:=-Y0;}
     {for i:=1 to nnt do
     begin
          smr:=smr+ro[i];
     end;
     r2:=smr/nnt;
     om:=-mu/r2;}
     {for i:=1 to nnt do
     begin
                x[i]:=(r2+0.1)*sin((om-0.001)*et[i]+om0+0.01)+x0+0.1;
                y[i]:=(r2+0.1)*cos((om-0.001)*et[i]+om0-0.01)+y0-0.1;
     end;}

           {assign(rs,'aflcirpr.res');
           rewrite(rs);
           for i:=1 to nnt do
           begin
                xc:=r2*sin(om*et[i]+om0)+x0;
                yc:=r2*cos(om*et[i]+om0)+y0;
                nevx:=x[i]-xc;
                nevy:=y[i]-yc;
                writeln(rs,(et[i]):8:3,x[i]:8:3,' ',y[i]:8:3,' ',nevx:8:3,' ',nevy:8:3,' ',ro[i]:8:3,' ',teta[i]:8:3);
                writeln(rs,(et[i]):8:3,xc:8:3,' ',yc:8:3,' ',nevx:8:3,' ',nevy:8:3,' ',ro[i]:8:3,' ',teta[i]:8:3);
           end;}
           {x0:=x0+r2*sin(om0)-r2*sin(om0+pi);
           y0:=y0+r2*cos(om0)-r2*cos(om0+pi);
           om0:=om0+pi;
           om:=-om;

           for i:=1 to nnt do
           begin
                xc:=r2*sin(om*et[i]+om0)+x0;
                yc:=r2*cos(om*et[i]+om0)+y0;
                nevx:=x[i]-xc;
                nevy:=y[i]-yc;
                writeln(rs,(et[i]):8:3,x[i]:8:3,' ',y[i]:8:3,' ',nevx:8:3,' ',nevy:8:3,' ',ro[i]:8:3,' ',teta[i]:8:3);
                writeln(rs,(et[i]):8:3,xc:8:3,' ',yc:8:3,' ',nevx:8:3,' ',nevy:8:3,' ',ro[i]:8:3,' ',teta[i]:8:3);
           end;}

           {for i:=1 to 1500 do
           begin
                xc:=r2*sin(om*i+om0)+x0;
                yc:=r2*cos(om*i+om0)+y0;
                writeln(rs,(i*1.0):8:3,xc:8:3,' ',yc:8:3);
           end;}

           {exit;}
     br2:=r2;
     bom:=om;
     bom0:=om0;
     bx0:=x0;
     by0:=y0;

     itr:=0;
     nsl:=false;

     l1:;
     repeat
           km:=0;
           smr:=0;
           for i:=1 to nnt do
           begin
                if gbop^[i]=0 then
                begin
                     km:=km+1;
                     sn:=sin(om*et^[i]+om0);
                     cs:=cos(om*et^[i]+om0);
                     a^[km  ,1]:=sn;
                     a^[km+1,1]:=cs;
                     a^[km  ,2]:= r2*cs*et^[i];
                     a^[km+1,2]:=-r2*sn*et^[i];
                     a^[km  ,3]:= r2*cs;
                     a^[km+1,3]:=-r2*sn;
                     a^[km  ,4]:= 1;
                     a^[km+1,4]:= 0;
                     a^[km  ,5]:= 0;
                     a^[km+1,5]:= 1;
                     xc:=r2*sn+x0;
                     yc:=r2*cs+y0;
                     nevx:=x^[i]-xc;
                     nevy:=y^[i]-yc;

                     a^[km  ,6]:= -nevx;
                     a^[km+1,6]:= -nevy;
                     km:=km+1;
                     smr:=smr+sqr(nevx)+sqr(nevy);
                 end;
           end;
           ptlst(a,xf,cv,km,5);
           {if itr=19 then
           begin
                itr:=19;
           end;}
           sgmnk(a,cv,xf,sgm,5,km);
           smn:=smr;
           if itr=359 then
           begin
                itr:=359;
           end;
           repeat
                  sml:=cirlam(l,r2,om,om0,x0,y0,xf[1],xf[2],xf[3],xf[4],xf[5]);
                  raz:=sml-smn;
                  if (raz>0) and (raz>1e-10) then l:=l*0.5;

           until not((raz>0) and (raz>1e-10)) or (l<1e-5);
           if l<0.1 then
           begin
                sml:=cirlam(1.2*l,r2,om,om0,x0,y0,xf[1],xf[2],xf[3],xf[4],xf[5]);
                if sml<smn then l:=l*1.1;
           end;
           smr:=smr/km;
           smr:=sqrt(smr);
           r2 := r2-l*xf[1];
           om := om-l*xf[2];
           om0:=om0-l*xf[3];
           x0 := x0-l*xf[4];
           y0 := y0-l*xf[5];
           vyh:=true;
           vyh:= vyh and (abs(xf[1])<1e-4);
           vyh:= vyh and (abs(xf[2])<1e-7);
           vyh:= vyh and (abs(xf[3])<1e-4);
           vyh:= vyh and (abs(xf[4])<1e-4);
           vyh:= vyh and (abs(xf[5])<1e-4);
           itr:=itr+1;
           {writeln(r2:10:5);}
           {if itr>1700 then
           begin
                vyh:=false;
           end;}
           {writeln(rs,itr:5,' ',r2:15:5,' ',om:25:15,' ',om0:25:15,' ',x0:15:5,' ',y0:15:5);
           flush(rs);}
           {if om*bom<0 then
           begin
                nsl:=true;
           end;}
           if (r2>3000) or (r2<0) then
           begin
                {if nsl then inud:=inud+1;}
                nsl:=true;
                r2:=br2;
                om:=bom;
                om0:=bom0;
                x0:=bx0;
                y0:=by0;

                {x0:=x0+r2*sin(om0)-r2*sin(om0+pi);
                y0:=y0+r2*cos(om0)-r2*cos(om0+pi);
                om0:=om0+pi;
                om:=-om;}


                l:=l/2;
                goto l1;
           end;

     until (vyh) or (itr>3000) {or (l<1e-5)};
     if vyh then
     begin
          write(rso,' 1',km:4,' ',r2:t13:t8,' ',om:t13:t8,' ',om0:t13:t8,' ',x0:t13:t8,' ',y0:t13:8,' ');
          write(rso,sgm[1]:t11:t8,' ',sgm[2]:t11:t8,' ',sgm[3]:t11:t8,' ',sgm[4]:t11:t8,' ',sgm[5]:t11:t8,' ');
          write(rso,sgm[6]:t11:t8);
     end else
     begin
          write(rso,' 0');
     end;
     flush(rso);
     if l<1e-5 then
     begin
          itr:=3002;
     end;
     if nsl then
     begin
          if itr>3000 then
          begin
               inud:=inud+1;
          end else
          begin
               iud:=iud+1;
          end;
     end;
     {close(rs);}

     itrf:=itr;
     if itr>3000 then
     begin
          vyh:=false;
     end;
           assign(rs,'flcirit.res');
           rewrite(rs);
           xc:=r2*sin(om0)+x0;
           yc:=r2*cos(om0)+y0;

           for i:=1 to nnt do
           begin
                xc:=r2*sin(om*et^[i]+om0)+x0;
                yc:=r2*cos(om*et^[i]+om0)+y0;
                nevx:=x^[i]-xc;
                nevy:=y^[i]-yc;
                writeln(rs,(et^[i]):8:3,x^[i]:8:3,' ',y^[i]:8:3,' ',nevx:8:3,' ',nevy:8:3,' ',ro^[i]:8:3,' ',teta^[i]:8:3);
                writeln(rs,(et^[i]):8:3,xc:8:3,' ',yc:8:3,' ',nevx:8:3,' ',nevy:8:3,' ',ro^[i]:8:3,' ',teta^[i]:8:3);
                writeln(rs,(et^[i]):8:3,xc-0.01:8:3,' ',yc-0.01:8:3,' ',nevx:8:3,' ',nevy:8:3,' ',ro^[i]:8:3,' ',teta^[i]:8:3);
           end;
           close(rs);



end;*)
procedure nachprib;
var r2:extended;
om,om0,x0,y0,om02:extended;
xa,xt,xtt,ya,yt,ytt,zn,r3:extended;
mom,mom0,mx0,my0,med,mr:extended;
tom,tom0,tx0,ty0,ted,nevx,nevy,epsr,x0s,y0s,x1s,y1s,x2s,y2s,ro1s,ro2s,teta1s,teta2s:extended;

begin
           (*assign(dt,'pvdt.res');
           rewrite(dt);
           for i:=1 to nnt do
           begin
                 writeln(dt,(et^[i]+t0):8:3,x^[i]:8:3,' ',y^[i]:8:3,' ',ro^[i]:8:3,' ',teta^[i]:8:3);
           end;
           flush(dt);
           close(dt);
           km:=0;
           for i:=t1 to nnt do
           begin
                km:=km+1;
                a^[km,1]:=et^[i]*et^[i];
                a^[km,2]:=et^[i];
                a^[km,3]:=1;
                a^[km,4]:=ro^[i];
           end;
           ptlst(a,xf,cv,km,3);
           sgmnk(a,cv,xf,sgm,3,km);
           svx:=sgm[2];
           nro0:=xf[3];
           sro0:=sgm[3];
           nrot:=xf[2];
           srot:=sgm[2];
           rot2:=xf[1];

           writeln('Oshibka ed. vesa po ro=',sgm[4]:8:4);
           writeln('Nevazki po ro=');
           km:=0;
           for i:=t1 to nnt do
           begin
                km:=km+1;
                writeln( et^[i]+t0:8:3,' ',(a^[km,4]-a^[km,3]*xf[3]-a^[km,2]*xf[2]-a^[km,1]*xf[1]):8:3);
           end;

           if brdln then  readln;

           km:=0;
           for i:=t1 to nnt do
           begin
                km:=km+1;
                a^[km,1]:=et^[i]*et^[i];
                a^[km,2]:=et^[i];
                a^[km,3]:=1;
                a^[km,4]:=teta^[i];
           end;
           ptlst(a,xf2,cv,km,3);
           sgmnk(a,cv,xf2,sgm,3,km);
           writeln('Nevazki po teta=');
           km:=0;
           for i:=t1 to nnt do
           begin
                km:=km+1;
                writeln('teta ', et^[i]+t0:8:3,' ',(a^[km,4]-a^[km,3]*xf2[3]-a^[km,2]*xf2[2]-a^[km,1]*xf2[1]):8:3);
           end;

           writeln('Oshibka ed. vesa po teta=',sgm[3]:8:4);
           ku:=ku+1;
           nteta0:=xf2[3];
           steta0:=sgm[3];
           ntetat:=xf2[2];

           stetat:=sgm[2];
           tetat2:=xf2[1];
           svy:=sgm[2];
           iteta0:=xf2[3];
           iro0:=xf[3];
           tetat:=xf2[2]{+xf2[1]/10};
           rot:=xf[2];
           xs:=iro0*sin(iteta0*pi/180);
           ys:=iro0*cos(iteta0*pi/180);
           vx:=rot*sin(iteta0*pi/180)+iro0*tetat*cos(iteta0*pi/180)*pi/180;
           vy:=rot*cos(iteta0*pi/180)-iro0*tetat*sin(iteta0*pi/180)*pi/180;
           mu:=sqrt(vx*vx+vy*vy);
           tetat:=tetat*pi/180;
           tetat2:=tetat2*pi/180;
           if bnkri then r:=abs(tetat*(mu*mu+rot*rot)+iro0*(rot*tetat2-rot2*tetat))/(mu*mu*mu);
           r2:=abs(tetat*(mu*mu+rot*rot)+iro0*(rot*tetat2-rot2*tetat))/(mu*mu*mu);
           r2:=1/r2;
           linprib(xa,xt,xtt,ya,yt,ytt,t1,nnt,a,cv,xf,sgm,x,y,et);
           linprib2(xa,xt,xtt,ya,yt,ytt,t1,nnt,a,cv,xf,sgm,x,y,et
           ,mr,mom,mom0,mx0,my0,posf);
           zn:=xt*xt+yt*yt;
           r3:=abs(xt*ytt-yt*xtt)/sqrt(zn*zn*zn);
           r3:=1/r3;
           {r3:=80;}
           dxw:=-vy;
           dyw:=vx;
           om0:=arctan(dxw/dyw);
           if dyw<0 then om0:=pi+om0;
           if om0<0 then om0:=om0+2*pi;
           {om0:=om0+pi/2;}
           {r2:=0;}
           med:=1e39;

           tom:=mu/r2;
           tom0:=om0;
           tx0:=xs-r2*sin(tom0);
           ty0:=ys-r2*cos(tom0);
           ted:=0;
           med:=1e39;
           for i:=1 to nnt do
           begin
                xc:=r2*sin(tom*et^[i]+tom0)+tx0;
                yc:=r2*cos(tom*et^[i]+tom0)+ty0;
                nevx:=x^[i]-xc;
                nevy:=y^[i]-yc;
                ted:=ted+nevx*nevx+nevy*nevy;
           end;
           ted:=sqrt(ted/(nnt*2));
           if med>ted then
           begin
                mr:=r2;
                med:=ted;
                mom:=tom;
                mom0:=tom0;
                mx0:=tx0;
                my0:=ty0;
           end;

           tom:=-mu/r2;
           tom0:=om0;
           tx0:=xs-r2*sin(tom0);
           ty0:=ys-r2*cos(tom0);
           ted:=0;
           for i:=1 to nnt do
           begin
                xc:=r2*sin(tom*et^[i]+tom0)+tx0;
                yc:=r2*cos(tom*et^[i]+tom0)+ty0;
                nevx:=x^[i]-xc;
                nevy:=y^[i]-yc;
                ted:=ted+nevx*nevx+nevy*nevy;
           end;
           ted:=sqrt(ted/(nnt*2));
           if med>ted then
           begin
                mr:=r2;
                med:=ted;
                mom:=tom;
                mom0:=tom0;
                mx0:=tx0;
                my0:=ty0;
           end;

           tom:=mu/r2;
           tom0:=om0+pi;
           tx0:=xs-r2*sin(tom0);
           ty0:=ys-r2*cos(tom0);
           ted:=0;
           for i:=1 to nnt do
           begin
                xc:=r2*sin(tom*et^[i]+tom0)+tx0;
                yc:=r2*cos(tom*et^[i]+tom0)+ty0;
                nevx:=x^[i]-xc;
                nevy:=y^[i]-yc;
                ted:=ted+nevx*nevx+nevy*nevy;
           end;
           ted:=sqrt(ted/(nnt*2));
           if med>ted then
           begin
                mr:=r2;
                med:=ted;
                mom:=tom;
                mom0:=tom0;
                mx0:=tx0;
                my0:=ty0;
           end;

           tom:=-mu/r2;
           tom0:=om0+pi;
           tx0:=xs-r2*sin(tom0);
           ty0:=ys-r2*cos(tom0);
           ted:=0;
           for i:=1 to nnt do
           begin
                xc:=r2*sin(tom*et^[i]+tom0)+tx0;
                yc:=r2*cos(tom*et^[i]+tom0)+ty0;
                nevx:=x^[i]-xc;
                nevy:=y^[i]-yc;
                ted:=ted+nevx*nevx+nevy*nevy;
           end;
           ted:=sqrt(ted/(nnt*2));
           if med>ted then
           begin
                mr:=r2;
                med:=ted;
                mom:=tom;
                mom0:=tom0;
                mx0:=tx0;
                my0:=ty0;
           end;

           dxw:=-yt;
           dyw:=xt;
           om0:=arctan(dxw/dyw);
           if dyw<0 then om0:=pi+om0;
           if om0<0 then om0:=om0+2*pi;
           mu:=sqrt(xt*xt+yt*yt);
           xs:=xa;
           ys:=ya;
           tom:=mu/r3;
           tom0:=om0;
           tx0:=xs-r3*sin(tom0);
           ty0:=ys-r3*cos(tom0);
           ted:=0;
           {med:=1e39;}
           for i:=1 to nnt do
           begin
                xc:=r3*sin(tom*et^[i]+tom0)+tx0;
                yc:=r3*cos(tom*et^[i]+tom0)+ty0;
                nevx:=x^[i]-xc;
                nevy:=y^[i]-yc;
                ted:=ted+nevx*nevx+nevy*nevy;
           end;
           ted:=sqrt(ted/(nnt*2));
           if med>ted then
           begin
                mr:=r3;
                med:=ted;
                mom:=tom;
                mom0:=tom0;
                mx0:=tx0;
                my0:=ty0;
           end;

           tom:=-mu/r3;
           tom0:=om0;
           tx0:=xs-r3*sin(tom0);
           ty0:=ys-r3*cos(tom0);
           ted:=0;
           for i:=1 to nnt do
           begin
                xc:=r3*sin(tom*et^[i]+tom0)+tx0;
                yc:=r3*cos(tom*et^[i]+tom0)+ty0;
                nevx:=x^[i]-xc;
                nevy:=y^[i]-yc;
                ted:=ted+nevx*nevx+nevy*nevy;
           end;
           ted:=sqrt(ted/(nnt*2));
           if med>ted then
           begin
                mr:=r3;
                med:=ted;
                mom:=tom;
                mom0:=tom0;
                mx0:=tx0;
                my0:=ty0;
           end;

           tom:=mu/r3;
           tom0:=om0+pi;
           tx0:=xs-r3*sin(tom0);
           ty0:=ys-r3*cos(tom0);
           ted:=0;
           for i:=1 to nnt do
           begin
                xc:=r3*sin(tom*et^[i]+tom0)+tx0;
                yc:=r3*cos(tom*et^[i]+tom0)+ty0;
                nevx:=x^[i]-xc;
                nevy:=y^[i]-yc;
                ted:=ted+nevx*nevx+nevy*nevy;
           end;
           ted:=sqrt(ted/(nnt*2));
           if med>ted then
           begin
                mr:=r3;
                med:=ted;
                mom:=tom;
                mom0:=tom0;
                mx0:=tx0;
                my0:=ty0;
           end;

           tom:=-mu/r3;
           tom0:=om0+pi;
           tx0:=xs-r3*sin(tom0);
           ty0:=ys-r3*cos(tom0);
           ted:=0;
           for i:=1 to nnt do
           begin
                xc:=r3*sin(tom*et^[i]+tom0)+tx0;
                yc:=r3*cos(tom*et^[i]+tom0)+ty0;
                nevx:=x^[i]-xc;
                nevy:=y^[i]-yc;
                ted:=ted+nevx*nevx+nevy*nevy;
           end;
           ted:=sqrt(ted/(nnt*2));
           if med>ted then
           begin
                mr:=r3;
                med:=ted;
                mom:=tom;
                mom0:=tom0;
                mx0:=tx0;
                my0:=ty0;
           end;*)
           linprib2(xa,xt,xtt,ya,yt,ytt,t1,nnt,a,cv,xf,sgm,x,y,et
           ,mr,mom,mom0,mx0,my0,posf,gbop,ae9);
           {flcir(mr,mom,mom0,mx0,my0);}


           epsr:=1e-6;
           x0s:=mr*sin(mom*0+mom0)+mx0;
           y0s:=mr*cos(mom*0+mom0)+my0;
           x1s:=mr*sin(mom*(-epsr)+mom0)+mx0;
           y1s:=mr*cos(mom*(-epsr)+mom0)+my0;
           x2s:=mr*sin(mom*(epsr)+mom0)+mx0;
           y2s:=mr*cos(mom*(epsr)+mom0)+my0;

           iro0:=sqrt(x0s*x0s+y0s*y0s);
           ro1s:=sqrt(x1s*x1s+y1s*y1s);
           ro2s:=sqrt(x2s*x2s+y2s*y2s);

           dxw:=x0s;
           dyw:=y0s;
           iteta0:=arctan(dxw/dyw);
           if dyw<0 then iteta0:=pi+iteta0;
           if iteta0<0 then iteta0:=iteta0+2*pi;
           iteta0:=iteta0*180/pi;

           dxw:=x1s;
           dyw:=y1s;
           teta1s:=arctan(dxw/dyw);
           if dyw<0 then teta1s:=pi+teta1s;
           if teta1s<0 then teta1s:=teta1s+2*pi;

           dxw:=x2s;
           dyw:=y2s;
           teta2s:=arctan(dxw/dyw);
           if dyw<0 then teta2s:=pi+teta2s;
           if teta2s<0 then teta2s:=teta2s+2*pi;
           rot:=0.5*(ro2s-ro1s)/epsr;
           tetat:=0.5*(teta2s-teta1s)/epsr;
           r:=mr;
           tetat:=tetat*180/pi;
           writeln('vx=',vx:10:5,'+/-',svx:8:5,' ','vy=',vy:10:5,'+/-',svy:8:5);
           writeln('vx=',vx:10:5,'+/-',svx:8:5,' ','vy=',vy:10:5,'+/-',svy:8:5);
           writeln('ro0=',iro0:10:5,' teta0=',iteta0:10:5,' tro:=',rot:10:5,' tetat:=',tetat:10:5);
          {r:=iro0;}
           writeln('tetat2=',tetat2:20:10);




          writeln('Krivizna=',r:8:3);
          if bnvr then
          begin
                vr:=sqrt((16/(pi*pi))-1)*mu;
                vr:=4.74*mu/pia;
          end;
          writeln('vr=',vr:10:5);
          writeln('massa=',massa:5:1);
          writeln('nt=',nt);
          if brdln then  readln;
          fd:=0;
          clrscr;

end;



(*procedure vyhtup(iro0,iteta0,rot,tetat,r,vr,massa,pia:extended);
const eps=1e-5;
      eps2=1e-6;

var
ti0,vgrad:das;
dro,dte,f0,f1,f2,l,itr,min,minl,f3,f00:extended;
jb,mi,itr2:integer;
begin
     f00:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
     itr2:=0;

     repeat
                     f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
                     jb:=0;
                     if bro0 then
                     begin
                          jb:=jb+1;
                          f1:=fnev(iro0+eps,iteta0,rot,tetat,r,vr,massa,pia);
                          vgrad[jb]:=(f1-f0)/eps;
                     end;

                     if bteta0 then
                     begin
                          dro:=eps*180/pi;
                          jb:=jb+1;
                          f1:=fnev(iro0,iteta0+dro,rot,tetat,r,vr,massa,pia);
                          vgrad[jb]:=(f1-f0)/dro;

                     end;
                     if brot then
                     begin
                          dro:=eps*0.01;
                          jb:=jb+1;
                          f1:=fnev(iro0,iteta0,rot+dro,tetat,r,vr,massa,pia);
                          vgrad[jb]:=(f1-f0)/dro;
                     end;
                     if btetat then
                     begin
                          jb:=jb+1;
                          f1:=fnev(iro0,iteta0,rot,tetat+eps,r,vr,massa,pia);
                          vgrad[jb]:=(f1-f0)/eps;

                     end;
                     if br then
                     begin
                          jb:=jb+1;
                          f1:=fnev(iro0,iteta0,rot,tetat,r+eps,vr,massa,pia);
                          vgrad[jb]:=(f1-f0)/eps;

                     end;


                     if bvr then
                     begin
                          jb:=jb+1;
                          f1:=fnev(iro0,iteta0,rot,tetat,r,vr+eps,massa,pia);
                          vgrad[jb]:=(f1-f0)/eps;

                     end;
                     if bmassa then
                     begin
                          jb:=jb+1;
                          f1:=fnev(iro0,iteta0,rot,tetat,r,vr,massa+eps,pia);
                          vgrad[jb]:=(f1-f0)/eps;

                     end;

                     if bpia then
                     begin
                          jb:=jb+1;
                          f1:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia+eps);
                          vgrad[jb]:=(f1-f0)/eps;
                     end;
                     jb:=0;
                     min:=1e39;
                     if bro0 then
                     begin
                          jb:=jb+1;
                          f2:=f0;
                          l:=sign(vgrad[jb])*1e-5;
                          itr:=0;
                          repeat
                                f1:=f2;
                                f2:=fnev(iro0-l,iteta0,rot,tetat,r,vr,massa,pia);
                                if f2<=f1 then l:=l*1.1;
                                itr:=itr+1;
                          until (f2>f1) or (itr>10000) ;
                          if itr>10000 then
                          begin
                               writeln('problema v procedure vyhtup');
                               halt(0);
                          end;
                          if f1<min then
                          begin
                               mi:=jb;
                               min:=f1;
                               minl:=l/1.1;
                          end;
                     end;

                     if bteta0 then
                     begin
                          jb:=jb+1;
                          f2:=f0;
                          l:=sign(vgrad[jb])*1e-5*180/pi;
                          itr:=0;
                          repeat
                                f1:=f2;
                                f2:=fnev(iro0,iteta0-l,rot,tetat,r,vr,massa,pia);
                                if f2<=f1 then l:=l*1.1;
                                itr:=itr+1;
                          until (f2>f1) or (itr>10000) ;
                          if itr>10000 then
                          begin
                               writeln('problema v procedure vyhtup');
                               halt(0);
                          end;
                          if f1<min then
                          begin
                               mi:=jb;
                               min:=f1;
                               minl:=l/1.1;
                          end;
                     end;
                     if brot then
                     begin
                          jb:=jb+1;

                          f2:=f0;
                          l:=sign(vgrad[jb])*0.01;
                          itr:=0;
                          repeat
                                f1:=f2;
                                f2:=fnev(iro0,iteta0,rot-l,tetat,r,vr,massa,pia);
                                if f2<=f1 then l:=l*1.1;
                                itr:=itr+1;
                          until (f2>f1) or (itr>10000) ;
                          if itr>10000 then
                          begin
                               writeln('problema v procedure vyhtup');
                               halt(0);
                          end;
                          if f1<min then
                          begin
                               mi:=jb;
                               min:=f1;
                               minl:=l/1.1;
                          end;
                     end;
                     if btetat then
                     begin
                          jb:=jb+1;
                          f2:=f0;
                          l:=sign(vgrad[jb]);
                          itr:=0;
                          repeat
                                f1:=f2;
                                f2:=fnev(iro0,iteta0,rot,tetat-l,r,vr,massa,pia);
                                if f2<=f1 then l:=l*1.1;
                                itr:=itr+1;
                          until (f2>f1) or (itr>10000) ;
                          if itr>10000 then
                          begin
                               writeln('problema v procedure vyhtup');
                               halt(0);
                          end;
                          if f1<min then
                          begin
                               mi:=jb;
                               min:=f1;
                               minl:=l/1.1;
                          end;
                     end;
                     if br then
                     begin
                          jb:=jb+1;
                          f2:=f0;
                          l:=sign(vgrad[jb])*1e-5;
                          itr:=0;
                          repeat
                                f1:=f2;
                                f2:=fnev(iro0,iteta0,rot,tetat,r-l,vr,massa,pia);
                                if f2<=f1 then l:=l*1.1;
                                itr:=itr+1;
                          until (f2>f1) or (itr>10000) ;
                          if itr>10000 then
                          begin
                               writeln('problema v procedure vyhtup');
                               halt(0);
                          end;
                          if f1<min then
                          begin
                               mi:=jb;
                               min:=f1;
                               minl:=l/1.1;
                          end;
                     end;


                     if bvr then
                     begin
                          jb:=jb+1;
                          f2:=f0;
                          l:=sign(vgrad[jb])*1e-5;
                          itr:=0;
                          repeat
                                f1:=f2;
                                f2:=fnev(iro0,iteta0,rot,tetat,r,vr-l,massa,pia);
                                if f2<=f1 then l:=l*1.1;
                                itr:=itr+1;
                          until (f2>f1) or (itr>10000) ;
                          if itr>10000 then
                          begin
                               writeln('problema v procedure vyhtup');
                               halt(0);
                          end;
                          if f1<min then
                          begin
                               mi:=jb;
                               min:=f1;
                               minl:=l/1.1;
                          end;

                     end;
                     if bmassa then
                     begin
                          jb:=jb+1;
                          f2:=f0;
                          l:=sign(vgrad[jb])*1e-5;
                          itr:=0;
                          repeat
                                f1:=f2;
                                f2:=fnev(iro0,iteta0,rot,tetat,r,vr,massa-l,pia);
                                if f2<=f1 then l:=l*1.1;
                                itr:=itr+1;
                          until (f2>f1) or (itr>10000) ;
                          if itr>10000 then
                          begin
                               writeln('problema v procedure vyhtup');
                               halt(0);
                          end;
                          if f1<min then
                          begin
                               mi:=jb;
                               min:=f1;
                               minl:=l/1.1;
                          end;
                     end;

                     if bpia then
                     begin
                          jb:=jb+1;
                          f2:=f0;
                          l:=sign(vgrad[jb])*1e-5;
                          itr:=0;
                          repeat
                                f1:=f2;
                                f2:=fnev(iro0,iteta0,rot,tetat-l,r,vr,massa,pia-l);
                                if f2<=f1 then l:=l*1.1;
                                itr:=itr+1;
                          until (f2>f1) or (itr>10000) ;
                          if itr>10000 then
                          begin
                               writeln('problema v procedure vyhtup');
                               halt(0);
                          end;
                          if f1<min then
                          begin
                               mi:=jb;
                               min:=f1;
                               minl:=l/1.1;
                          end;
                     end;
                     jb:=0;
                     if bro0 then
                     begin
                          jb:=jb+1;
                          if mi=jb then iro0:=iro0-minl;
                     end;

                     if bteta0 then
                     begin
                          jb:=jb+1;
                          if mi=jb then iteta0:=iteta0-minl;

                     end;
                     if brot then
                     begin
                          jb:=jb+1;
                          if mi=jb then rot:=rot-minl;

                     end;
                     if btetat then
                     begin
                          jb:=jb+1;
                          if mi=jb then tetat:=tetat-minl;
                     end;
                     if br then
                     begin
                          jb:=jb+1;
                          if mi=jb then r:=r-minl;
                     end;
                     if bvr then
                     begin
                          jb:=jb+1;
                          if mi=jb then vr:=vr-minl;

                     end;
                     if bmassa then
                     begin
                          jb:=jb+1;
                          if mi=jb then massa:=massa-minl;

                     end;

                     if bpia then
                     begin
                          jb:=jb+1;
                          if mi=jb then pia:=pia-minl;
                     end;
                     f3:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);;
                 itr2:=itr2+1;
           until itr2>100;


           f2:=0;


end;*)

procedure maxnev2(iro0,iteta0,rot,tetat,r,vr,massa,pia:extended;var oevr,oevt:extended);
var
ti0:das;
dro,dte,sum2,sum:extended;
begin
      {nt:=1;}
      mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);
      {if ti0[2]>=1 then
      begin
           fnev:=1e39;
      end else}
      begin

      sum:=0;
      sum2:=0;
      for id:=1 to nt do
      begin
           mpdq.rtti(ti0,et^[id],rog,tetag);
           dro:=rog-ro^[id];
           dte:=(tetag-teta^[id]){*pi*ro[id]/180};
           ogr2(dte);
           sum:=sum+dro*dro;
           sum2:=sum2+dte*dte;
           {sum:=sum+abs(dro)+abs(dte);}
           {sum:=sum+sqrt(abs(dro))+sqrt(abs(dte));}


      end;
      sum:=sum/nt;
      sum:=sqrt(sum);
      sum2:=sum2/nt;
      sum2:=sqrt(sum2);
      oevr:=sum;
      oevt:=sum2;
      end;
end;


(*procedure iskl2(iro0,iteta0,rot,tetat,r,vr,massa,pia,oevr,oevt:extended;var maxnev:extended;var mi:integer);
var
ti0:das;
dro,dte:extended;
begin
      {nt:=1;}
      mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);
      {if ti0[2]>=1 then
      begin
           fnev:=1e39;
      end else}
      begin

      maxnev:=0;
      mi:=-1;
      for id:=1 to nt do
      begin
           mpdq.rtti(ti0,et^[id],rog,tetag);
           dro:=abs(rog-ro^[id])/oevr;
           dte:=abs(tetag-teta^[id]){/oevt};
           ogr2(dte);
           dte:=abs(dte)/oevt;
           if dte>dro then dro:=dte;
           if dro>maxnev then
           begin
                maxnev:=dro;
                mi:=id;
           end;
           roo^[id]:=rog;
           tetao^[id]:=tetag;



      end;
      end;
end;*)
function fnev2(iro0,iteta0,rot,tetat,r,vr,massa,pia:extended):extended;
var
ti0:das;
dro,dte:extended;
begin
      mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);
      ti0[2]:=ti0[2]+0.0001;
      {if ti0[2]>=1 then
      begin
           fnev:=1e39;
      end else}
      begin

      sum:=0;
      for id:=1 to nt do
      begin
           mpdq.rtti(ti0,et^[id],rog,tetag);
           dro:=rog-ro^[id];
           dte:=(tetag-teta^[id])*pi*ro^[id]/180;
           sum:=sum+dro*dro+dte*dte;
           {sum:=sum+abs(dro)+abs(dte);}
           {sum:=sum+sqrt(abs(dro))+sqrt(abs(dte));}


      end;
      sum:=0.5*sum/nt;
      sum:=sqrt(sum);
      fnev2:=sum;
      end;
end;

(*procedure iteracii3(biskl:boolean);
const eps=1e-5;
      eps2=1e-6;

var itr,itr2:longint;
    vgrad:das;
    f0,f1,y,lam,f2,f3,prf0,f4,f5,l,l0,eps3,fs,f6,lam2,lam3:extended;
    ch:char;

BEGIN
      dro:=eps*180/pi;
      f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
      f1:=fnev(iro0+eps,iteta0+dro,rot+eps,tetat+eps,r+eps,vr+eps,massa,pia);

      itr:=0;
      lam:=1e-8;
      gllam:=1e-4;
      prf0:=1e6;

      repeat
             itr:=itr+1;
             f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
             jb:=0;
                     if bro0 then
                     begin
                          jb:=jb+1;
                          f1:=fnev(iro0+eps,iteta0,rot,tetat,r,vr,massa,pia);
                          vgrad[jb]:=(f1-f0)/eps;
                     end;

                     if bteta0 then
                     begin
                          dro:=eps*180/pi;
                          jb:=jb+1;
                          f1:=fnev(iro0,iteta0+dro,rot,tetat,r,vr,massa,pia);
                          vgrad[jb]:=(f1-f0)/dro;

                     end;
                     if brot then
                     begin
                          dro:=eps*0.01;
                          jb:=jb+1;
                          f1:=fnev(iro0,iteta0,rot+dro,tetat,r,vr,massa,pia);
                          vgrad[jb]:=(f1-f0)/dro;
                     end;
                     if btetat then
                     begin
                          jb:=jb+1;
                          f1:=fnev(iro0,iteta0,rot,tetat+eps,r,vr,massa,pia);
                          vgrad[jb]:=(f1-f0)/eps;

                     end;
                     if br then
                     begin
                          jb:=jb+1;
                          f1:=fnev(iro0,iteta0,rot,tetat,r+eps,vr,massa,pia);
                          vgrad[jb]:=(f1-f0)/eps;

                     end;


                     if bvr then
                     begin
                          jb:=jb+1;
                          f1:=fnev(iro0,iteta0,rot,tetat,r,vr+eps,massa,pia);
                          vgrad[jb]:=(f1-f0)/eps;

                     end;
                     if bmassa then
                     begin
                          jb:=jb+1;
                          f1:=fnev(iro0,iteta0,rot,tetat,r,vr,massa+eps,pia);
                          vgrad[jb]:=(f1-f0)/eps;

                     end;

                     if bpia then
                     begin
                          jb:=jb+1;
                          f1:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia+eps);
                          vgrad[jb]:=(f1-f0)/eps;
                     end;
                 lam:=abs(lam);
                f2:=flam(lam,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);

                while f2>f0 do
                begin
                     lam:=lam/1.1;
                     f2:=flam(lam,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                end;
                if lam<1e-10 then
                begin
                     lam3:=-1e-15;
                     f3:=f0;
                     repeat
                          f2:=f3;
                          f3:=flam(lam3,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                          lam3:=lam3*1.1;
                     until f2<f3;
                     if f2<f0 then
                     begin
                          lam:=lam3/1.1;
                     end;
                end;
                {while f2>=f0 do
                begin
                     lam3:=lam3/1.1;
                     f2:=flam(lam3,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                end;}


             f3:=flam(lam/1.1,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
             f4:=flam(lam*1.03,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
             if f3<f2 then lam:=lam/1.1;
             if f4<f2 then lam:=lam*1.03;
             l:=0;
             {f2:=flam(lam,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
             f3:=fff(lam,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
             l:=0;
             itr2:=0;
             lam2:=0.1;
             repeat

                   f4:=fff(l,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                   eps3:=gllam*0.0001;
                   f6:=fff(l-eps3*0.5,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);

                   f5:=fff(l+eps3*0.5,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                   fs:=(f5-f6)/eps3;
                   l0:=l;
                   l:=l-(lam2*f4/fs);
                   itr2:=itr2+1;
             until abs(f4)<1;
             l:=l*0.5;
             lam:=l;
             gllam:=l;}

                     jb:=0;
                     if bro0 then
                     begin
                          jb:=jb+1;
                          iro0:=iro0-lam*vgrad[jb];
                     end;

                     if bteta0 then
                     begin
                          jb:=jb+1;
                          iteta0:=iteta0-lam*vgrad[jb];
                     end;
                     if brot then
                     begin
                          jb:=jb+1;
                          rot:=rot-lam*vgrad[jb];

                     end;
                     if btetat then
                     begin
                          jb:=jb+1;
                          tetat:=tetat-lam*vgrad[jb];
                     end;
                     if br then
                     begin
                          jb:=jb+1;
                          r:=r-lam*vgrad[jb];
                     end;


                     if bvr then
                     begin
                          jb:=jb+1;
                          vr:=vr-lam*vgrad[jb];
                     end;
                     if bmassa then
                     begin
                          jb:=jb+1;
                          massa:=massa-lam*vgrad[jb];
                     end;

                     if bpia then
                     begin
                          jb:=jb+1;
                          pia:=pia-lam*vgrad[jb];
                     end;
                     y:=0;
                if itr mod 1 =0 then
                begin
                gotoxy(1,1);
                writeln;
                if bro0 then pst:='     var'
                        else pst:=' not var'  ;
                writeln('   ro0=',iro0:10:5,pst);
                if bteta0 then pst:='     var'
                        else pst:=' not var'  ;
                writeln(' teta0=',iteta0:10:5,pst);
                if brot then pst:='     var'
                        else pst:=' not var'  ;
                writeln('   rot=',rot:20:10,pst);
                if btetat then pst:='     var'
                        else pst:=' not var'  ;
                writeln(' tetat=',tetat:10:5,pst);
                if br then pst:='     var'
                        else pst:=' not var'  ;
                writeln('krivizn=',r:10:5,pst);
                if bvr then pst:='     var'
                        else pst:=' not var'  ;
                writeln('     vr=',vr:10:5,pst);
                if bmassa then pst:='     var'
                         else pst:=' not var'  ;
                writeln('  massa=',massa:10:5,pst);
                if bpia then pst:='     var'
                         else pst:=' not var'  ;
                writeln('     pi=',pia:10:5,pst);


                writeln('iteracija=',itr,' chislo per=',jb);

                writeln('f0(oev)=',f0:10:5);
                writeln('lam=',lam:20:15);
                writeln('e=',ev:20:15);
                wvh;
                end;
                y:=0;
            {if (f0>=prf0) then lam:=lam*0.5;}
            prf0:=f0;
            ch:='p';
            if keypressed then ch:=readkey;


      until (ch='a') or (itr>10000000);
END;*)
procedure ptir(ti0,tiro0,tiro,titeta0,titeta,
               tirot,tirot0,titetat0,titetat,
               tir0,tir,tivr0,tivr,tim0,tim,tip0,tip
               :das;
               var rog0,tetag0,eps,raz:extended;
               var i:integer;
               var jb,id:integer;
               et,ro,teta:dst3;
               var rog,tetag,rog1,tetag1,dro:extended;
               bro0,bteta0,brot,btetat,br,bvr,bmassa,bpia:boolean;
               a:msls);
begin

                     i:=i+1;
                     jb:=0;
                     mpdq.rtti(ti0,et^[id],rog,tetag);
                     {if id=1 then
                     begin
                          writeln(otl,'iro0=',iro0:10:5);
                          writeln(otl,'r=',r:10:5);
                          writeln(otl,'teta0=',iteta0:10:5);
                          writeln(otl,'rot=',rot:10:5);
                          writeln(otl,'tetat=',tetat:10:5);
                          writeln(otl,'massa=',massa:10:5);
                          writeln(otl,'vr=',vr:10:5);
                    end;}
                     if bro0 then
                     begin
                          mpdq.rtti(tiro0,et^[id],rog0,tetag0);
                          mpdq.rtti(tiro,et^[id],rog1,tetag1);
                          a^[i*2-1,jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[i]/(eps*180);
                          jb:=jb+1;
                     end;

                     if bteta0 then
                     begin
                          mpdq.rtti(titeta0,et^[id],rog0,tetag0);
                          mpdq.rtti(titeta,et^[id],rog1,tetag1);
                          dro:=eps*180/pi;
                          a^[i*2-1,jb+1]:=(rog1-rog0)/dro;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[i]/(dro*180);
                          jb:=jb+1;
                     end;
                     if brot then
                     begin
                          mpdq.rtti(tirot0,et^[id],rog0,tetag0);
                          mpdq.rtti(tirot,et^[id],rog1,tetag1);
                          dro:=eps*0.01;
                          a^[i*2-1,jb+1]:=(rog1-rog0)/dro;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[i]/(dro*180);
                          jb:=jb+1;
                     end;
                     if btetat then
                     begin
                          mpdq.rtti(titetat0,et^[id],rog0,tetag0);
                          mpdq.rtti(titetat,et^[id],rog1,tetag1);
                          a^[i*2-1,jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[i]/(eps*180);
                          jb:=jb+1;
                     end;
                     if br then
                     begin
                          mpdq.rtti(tir0,et^[id],rog0,tetag0);
                          mpdq.rtti(tir,et^[id],rog1,tetag1);
                          a^[i*2-1,jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[i]/(eps*180);
                          jb:=jb+1;
                     end;


                     if bvr then
                     begin
                          mpdq.rtti(tivr0,et^[id],rog0,tetag0);
                          mpdq.rtti(tivr,et^[id],rog1,tetag1);
                          a^[i*2-1,jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[i]/(eps*180);
                          jb:=jb+1;
                     end;
                     if bmassa then
                     begin
                          mpdq.rtti(tim0,et^[id],rog0,tetag0);
                          mpdq.rtti(tim,et^[id],rog1,tetag1);
                          a^[i*2-1,jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[i]/(eps*180);
                          jb:=jb+1;
                     end;

                     if bpia then
                     begin
                          mpdq.rtti(tip0,et^[id],rog0,tetag0);
                          mpdq.rtti(tip,et^[id],rog1,tetag1);
                          a^[i*2-1,jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[i]/(eps*180);

                          jb:=jb+1;
                     end;
                     a^[i*2-1,jb+1]:=rog-ro^[i];
                     raz:=(tetag-teta^[id]);
                     ogr2(raz);
                     a^[i*2,jb+1]:=raz*pi*ro^[i]/180;
end;



(*procedure iteracii(biskl:boolean);
var i,ih:integer;
    lp,rog2,ad,bd,h,rog0,tetag0:extended;
    rsx:text;
    ti0,tiro,titeta,tirot,titetat,tir,tivr,tim,tip,vgrad:das;
    tiro0,titeta0,tirot0,titetat0,tir0,tivr0,tim0,tip0:das;

    ch:char;
    f0,f1,f2,f3,f4,maxnev,ff1,ff2,raz:extended;
    mi,itr:integer;
    eps,eps2:extended;




begin
      eps:=1e-7;
      eps2:=1e-6;
      l:=1e-5;



      dro:=eps*180/pi;
      {f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
      f1:=fnev(iro0,iteta0,rot+eps,tetat,r,vr,massa,pia);}

      {f2:=fnev(iro0,iteta0,rot,tetat,r,vr-0.001,massa,pia);
      f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
      f3:=fnev(iro0,iteta0,rot,tetat,r,vr+0.001,massa,pia);}


      ch:='p';

      {h:=fnev;}
         fd:=0;
                {str(fd mod 2,pst);}
                {assign(rsx,'ret.res');
                rewrite(rsx);
                ad:=0.5539;
                bd:=2.5541;
                h:=(bd-ad)/10000;

                for id:=0 to 10000 do
                begin
                     dro:=ad+h*id;
                     mpdq.tpvd(iro0,iteta0,rot,tetat,dro,vr,massa,et[1],pia,ti0);
                     mpdq.rtti(ti0,et[1],rog,tetag);

                     writeln(rsx,id,' ',dro,' ',rog,' ',tetag);
                end;
                close(rsx);}

          begin
                {massa:=massa-0.1;}
          repeat
          t3:=nt;


          repeat
                i:=0;
                mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);
                ev:=ti0[2];
                {if ti0[2]>1 then
                begin
                      exit;
                end;}
                {if r<0.0001 then exit;}
                {if fd>550 then exit;}
        if bro0 then mpdq.tpvd(iro0-0.5*eps,iteta0,rot,tetat,r,vr,massa,et^[1],pia,tiro0);
        if bteta0 then mpdq.tpvd(iro0,iteta0-0.5*eps*180/pi,rot,tetat,r,vr,massa,et^[1],pia,titeta0);
        if brot then mpdq.tpvd(iro0,iteta0,rot-0.5*eps*0.01,tetat,r,vr,massa,et^[1],pia,tirot0);
        if btetat then mpdq.tpvd(iro0,iteta0,rot,tetat-0.5*eps,r,vr,massa,et^[1],pia,titetat0);
        if br then mpdq.tpvd(iro0,iteta0,rot,tetat,r-0.5*eps,vr,massa,et^[1],pia,tir0);
        if bvr then mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr-0.5*eps,massa,et^[1],pia,tivr0);
        if bmassa then mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa-0.5*eps,et^[1],pia,tim0);
        if bpia then  mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia-0.5*eps,tip0);

        if bro0 then mpdq.tpvd(iro0+0.5*eps,iteta0,rot,tetat,r,vr,massa,et^[1],pia,tiro);
        if bteta0 then mpdq.tpvd(iro0,iteta0+0.5*eps*180/pi,rot,tetat,r,vr,massa,et^[1],pia,titeta);
        if brot then mpdq.tpvd(iro0,iteta0,rot+0.5*eps*0.01,tetat,r,vr,massa,et^[1],pia,tirot);
        if btetat then mpdq.tpvd(iro0,iteta0,rot,tetat+0.5*eps,r,vr,massa,et^[1],pia,titetat);
        if br then mpdq.tpvd(iro0,iteta0,rot,tetat,r+0.5*eps,vr,massa,et^[1],pia,tir);
        if bvr then mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr+0.5*eps,massa,et^[1],pia,tivr);
        if bmassa then mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa+0.5*eps,et^[1],pia,tim);
        if bpia then  mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia+0.5*eps,tip);

                i:=0;
                for id:=1 to t3 do
                begin
                     if gbop^[id]=0 then
                     begin
                          ptir(ti0,tiro0,tiro,titeta0,titeta,
                               tirot,tirot0,titetat0,titetat,
                                tir0,tir,tivr0,tivr,tim0,tim,tip0,tip,
                                rog0,tetag0,eps,raz,i,
                                jb,id,
                                et,ro,teta,
                                rog,tetag,rog1,tetag1,dro,
                                bro0,bteta0,brot,btetat,br,bvr,bmassa,bpia,
                                a
                                 );
                     end;
                end;
                {halt(0);}
                {close(rsx);}
                ptlst(a,xf,cv,i*2,jb);
                {sgmnk(a,cv,xf,sgm,jb,i*2);}

                for ih:=1 to jb do vgrad[ih]:=xf[ih];
                f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);

                f2:=flam(l,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                l:=abs(l);
                while f2>f0 do
                begin
                     l:=l/1.1;

                     f2:=flam(l,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                end;

                {if abs(l)<1e-8 then
                begin
                     l:=-l;
                          ff2:=f0;
                          itr:=0;
                          repeat
                                ff1:=ff2;
                                ff2:=flam(l,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                                if ff2<=ff1 then l:=l*1.1;
                                itr:=itr+1;
                          until (ff2>ff1) or (itr>10000) ;
                          l:=l/1.1;

                end else}
                begin
                     f3:=flam(l/1.1,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                     f4:=flam(l*1.03,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                     if f3<f2 then l:=l/1.1;
                     if f4<=f2 then l:=l*1.03;
                end;


                     jb:=0;
                     vyh:=true;

                     lp:=l;
                     if bro0 then
                     begin
                          jb:=jb+1;
                          iro0:=iro0-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<eps2);
                     end;

                     if bteta0 then
                     begin
                          jb:=jb+1;
                          iteta0:=iteta0-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<eps2);
                     end;
                     if brot then
                     begin

                          jb:=jb+1;
                          rot:=rot-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<eps2);
                     end;
                     if btetat then
                     begin
                          jb:=jb+1;
                          tetat:=tetat-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb]/tetat)<eps2);

                     end;
                     if br then
                     begin
                          jb:=jb+1;
                          r:=r-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb]/r)<eps2);

                     end;

                     if bvr then
                     begin
                          jb:=jb+1;;
                          vr:=vr-lp*xf[jb];

                          vyh:=vyh and (abs(xf[jb]/vr)<eps2);

                     end;

                     if bmassa then
                     begin
                          jb:=jb+1;
                          massa:=massa-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb]/massa)<eps2);

                     end;

                     if bpia then
                     begin
                          jb:=jb+1;
                          pia:=pia-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb]/pia)<eps2);
                     end;
                fd:=fd+1;
                {if r<0.1 then exit;
                if r>15 then exit;}
                if massa<0 then exit;



                if fd mod 1 =0 then
                begin
                gotoxy(1,1);
                writeln;
                if bro0 then pst:='     var'
                        else pst:=' not var'  ;
                writeln('   ro0=',iro0:10:5,pst);
                if bteta0 then pst:='     var'
                        else pst:=' not var'  ;
                writeln(' teta0=',iteta0:10:5,pst);
                if brot then pst:='     var'
                        else pst:=' not var'  ;
                writeln('   rot=',rot:20:10,pst);
                if btetat then pst:='     var'
                        else pst:=' not var'  ;
                writeln(' tetat=',tetat:10:5,pst);
                if br then pst:='     var'
                        else pst:=' not var'  ;
                writeln('krivizn=',r:10:5,pst);
                if bvr then pst:='     var'
                        else pst:=' not var'  ;
                writeln('     vr=',vr:10:5,pst);
                if bmassa then pst:='     var'
                         else pst:=' not var'  ;
                writeln('  massa=',massa:10:5,pst);
                if bpia then pst:='     var'
                         else pst:=' not var'  ;
                writeln('     pi=',pia:10:5,pst);


                writeln('iteracija=',fd,' kol-vo urav=',i,' iz=',iz );

                writeln('f0(oev)=',f0:10:5 );
                writeln('lp=',lp );
                writeln('e=',ev:10:5 );

                if biskl then wvh;
                {readln;}
                if fd=172 then
                begin
                     fd:=172;
                end;

                {if fd mod 10 =0 then
                begin
                     ch:=readkey;
                     if ch='q' then halt(0);
                end;}
                {vydres(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt);;}

                end;
                {writeln(rsq,fd:5,' ',f0:10:5,' ',massa:10:5);
                flush(rsq);}

                {if fd mod 100 =0 then if biskl then vydres(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,
                nt,sgm,bro0,bteta0,brot,btetat,br
                ,bvr,bmassa,bpia
                ,nro0,sro0,nteta0,steta0,nrot,srot,ntetat,stetat,oevr,oevt);;}
             ch:='p';
            if keypressed then ch:=readkey;

            if (l<1e-10) or (ev>50) then
            begin
                 {vyhtup(iro0,iteta0,rot,tetat,r,vr,massa,pia);}
                 vyh:=true;
            end;



          until (ch='a') or vyh;
          {until vyh;}
          {vydres(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt);;;;}

          {nevotap(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt);;;;}
          if (l>1e-10) and (ev<50) then
          begin

          if biskl then maxnev2(iro0,iteta0,rot,tetat,r,vr,massa,pia,oevr,oevt);
          if biskl then
          begin
               iskl2(iro0,iteta0,rot,tetat,r,vr,massa,pia,oevr,oevt,maxnev,mi);
               ew:=maxnev>porog;
               if ew then
               begin
                    cvb(mi);
               end;


          end
          else ew:=false;
          {ew:=false;}
          {if biskl then wvh;}
          end;
          wvh;
          if (l<1e-10) or (ev>50) then ch:='a';


          until  (ch='a') or not(ew);
          {if biskl then vydres2(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt,oevr,oevt);}
          end;

end;*)
{procedure tcir(ti0,tiro0,tiro,titeta0,titeta,
               tirot,tirot0,titetat0,titetat,
               tir0,tir,tivr0,tivr,tim0,tim,tip0,tip
               :das;
               var rog0,tetag0,eps,raz:extended;
               var i:integer);
begin
                     jb:=0;
                     i:=i+1;
                     mpdq.rtti(ti0,et^[id],rog,tetag);
                     if bro0 then
                     begin
                          mpdq.rtti(tiro0,et^[id],rog0,tetag0);
                          mpdq.rtti(tiro,et^[id],rog1,tetag1);
                          a^[i*2-1,jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;

                     if bteta0 then
                     begin
                          mpdq.rtti(titeta0,et^[id],rog0,tetag0);
                          mpdq.rtti(titeta,et^[id],rog1,tetag1);
                          dro:=eps*180/pi;
                          a^[i*2-1,jb+1]:=(rog1-rog0)/dro;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[id]/(dro*180);
                          jb:=jb+1;
                     end;
                     if brot then
                     begin
                          mpdq.rtti(tirot0,et^[id],rog0,tetag0);
                          mpdq.rtti(tirot,et^[id],rog1,tetag1);
                          dro:=eps*0.01;
                          a^[i*2-1,jb+1]:=(rog1-rog0)/dro;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[id]/(dro*180);
                          jb:=jb+1;
                     end;
                     if btetat then
                     begin
                          mpdq.rtti(titetat0,et^[id],rog0,tetag0);
                          mpdq.rtti(titetat,et^[id],rog1,tetag1);
                          a^[i*2-1,jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;
                     if br then
                     begin
                          mpdq.rtti(tir0,et^[id],rog0,tetag0);
                          mpdq.rtti(tir,et^[id],rog1,tetag1);
                          a^[i*2-1,jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;


                     if bvr then
                     begin
                          mpdq.rtti(tivr0,et^[id],rog0,tetag0);
                          mpdq.rtti(tivr,et^[id],rog1,tetag1);
                          a^[i*2-1,jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;
                     if bmassa then
                     begin
                          mpdq.rtti(tim0,et^[id],rog0,tetag0);
                          mpdq.rtti(tim,et^[id],rog1,tetag1);
                          a^[i*2-1,jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;

                     if bpia then
                     begin
                          mpdq.rtti(tip0,et^[id],rog0,tetag0);
                          mpdq.rtti(tip,et^[id],rog1,tetag1);
                          a^[i*2-1,jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2,jb+1]:=raz*pi*ro^[id]/(eps*180);

                          jb:=jb+1;
                     end;
                     a^[i*2-1,jb+1]:=rog-ro^[id];
                     raz:=(tetag-teta^[id]);
                     ogr2(raz);
                     a^[i*2,jb+1]:=raz*pi*ro^[id]/180;

end;}


procedure proizv;


const eps=1e-5;
      eps2=1e-6;
var i:integer;
    h,t:extended;
    ti0,tiro,titeta,tirot,titetat,tir,tivr,tim,tip:das;
     resp:text;

begin
     mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);

     mpdq.tpvd(iro0+eps,iteta0,rot,tetat,r,vr,massa,et^[1],pia,tiro);
     mpdq.tpvd(iro0,iteta0+0.01*eps*180/pi,rot,tetat,r,vr,massa,et^[1],pia,titeta);
     mpdq.tpvd(iro0,iteta0,rot+eps*0.01,tetat,r,vr,massa,et^[1],pia,tirot);
     mpdq.tpvd(iro0,iteta0,rot,tetat+eps,r,vr,massa,et^[1],pia,titetat);
     mpdq.tpvd(iro0,iteta0,rot,tetat,r+eps,vr,massa,et^[1],pia,tir);
     mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr+eps,massa,et^[1],pia,tivr);
     mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa+eps,et^[1],pia,tim);
     mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia+eps,tip);
     h:=(et^[nt]-et^[1])/1000;

     assign(resp,'przvro.res');
     rewrite(resp);
     for i:=0 to 1000 do
     begin
          t:=et^[1]+i*h;
          mpdq.rtti(ti0,t,rog,tetag);
          mpdq.rtti(tiro,t,rog1,tetag1);
          dro:=(rog1-rog)/eps;
          dteta:=(tetag1-tetag)*pi*rog/(eps*180);
          writeln(resp,t+t0,' ',dro,' ',dteta);
     end;
     close(resp);

     assign(resp,'przvteta.res');
     rewrite(resp);
     for i:=0 to 1000 do
     begin
          t:=et^[1]+i*h;
          mpdq.rtti(ti0,t,rog,tetag);
          mpdq.rtti(titeta,t,rog1,tetag1);
          dro:=(rog1-rog)/(eps*0.01);
          dteta:=(tetag1-tetag)*pi*rog/(eps*180*0.01);
          writeln(resp,t+t0,' ',dro,' ',dteta);
     end;
     close(resp);

     assign(resp,'przvrot.res');
     rewrite(resp);
     for i:=0 to 1000 do
     begin
          t:=et^[1]+i*h;
          mpdq.rtti(ti0,t,rog,tetag);
          mpdq.rtti(tirot,t,rog1,tetag1);
          dro:=(rog1-rog)/eps;
          dteta:=(tetag1-tetag)*pi*rog/(eps*180);
          writeln(resp,t+t0,' ',dro,' ',dteta);
     end;
     close(resp);

     assign(resp,'przvtt.res');
     rewrite(resp);
     for i:=0 to 1000 do
     begin
          t:=et^[1]+i*h;
          mpdq.rtti(ti0,t,rog,tetag);
          mpdq.rtti(titetat,t,rog1,tetag1);
          dro:=(rog1-rog)/eps;
          dteta:=(tetag1-tetag)*pi*rog/(eps*180);
          writeln(resp,t+t0,' ',dro,' ',dteta);
     end;
     close(resp);

     assign(resp,'przvrc.res');
     rewrite(resp);
     for i:=0 to 1000 do
     begin
          t:=et^[1]+i*h;
          mpdq.rtti(ti0,t,rog,tetag);
          mpdq.rtti(tir,t,rog1,tetag1);
          dro:=(rog1-rog)/eps;
          dteta:=(tetag1-tetag)*pi*rog/(eps*180);
          writeln(resp,t+t0,' ',dro,' ',dteta);
     end;
     close(resp);

     assign(resp,'przvvr.res');
     rewrite(resp);
     for i:=0 to 1000 do
     begin
          t:=et^[1]+i*h;
          mpdq.rtti(ti0,t,rog,tetag);
          mpdq.rtti(tivr,t,rog1,tetag1);
          dro:=(rog1-rog)/eps;
          dteta:=(tetag1-tetag)*pi*rog/(eps*180);
          writeln(resp,t+t0,' ',dro,' ',dteta);
     end;
     close(resp);

     assign(resp,'przvm.res');
     rewrite(resp);
     for i:=0 to 1000 do
     begin
          t:=et^[1]+i*h;
          mpdq.rtti(ti0,t,rog,tetag);
          mpdq.rtti(tim,t,rog1,tetag1);
          dro:=(rog1-rog)/eps;
          dteta:=(tetag1-tetag)*pi*rog/(eps*180);
          writeln(resp,t+t0,' ',dro,' ',dteta);
     end;
     close(resp);

     assign(resp,'przvp.res');
     rewrite(resp);
     for i:=0 to 1000 do
     begin
          t:=et^[1]+i*h;
          mpdq.rtti(ti0,t,rog,tetag);
          mpdq.rtti(tip,t,rog1,tetag1);
          dro:=(rog1-rog)/eps;
          dteta:=(tetag1-tetag)*pi*rog/(eps*180);
          writeln(resp,t+t0,' ',dro,' ',dteta);
     end;
     close(resp);
end;

(*procedure iteracii2(biskl:boolean);
var i:integer;
    lp,rog2,x0,y0,xp,yp:extended;
    rsx:text;
    xf2,xf1,xf0:stb;
    kvyh:boolean;
begin
          begin
                {massa:=massa-0.1;}
          repeat
          t3:=nt;
          lp:=l;
          for i:=1 to 7 do xf1[i]:=0;
          for i:=1 to 7 do xf2[i]:=0;
          repeat
                i:=0;
                {str(fd mod 2,pst);
                assign(rsx,'ret'+pst+'.res');
                rewrite(rsx);
                for id:=-5000 to 5000 do
                begin
                     dro:=id*0.0001;
                     mpdq.tpvd(iro0,iteta0,rot,tetat,r+dro,vr,massa,et[300],pia,rog1,tetag1,t0,ev);
                     writeln(rsx,id,' ',r+dro,' ',rog1,' ',tetag1);
                end;}


                for id:=1 to t3 do
                begin
                     jb:=0;
                     i:=id;
                     mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et[id],pia,rog,tetag,t0,ev);
                     x0:=rog*sin(tetag*pi/180);
                     y0:=rog*cos(tetag*pi/180);

                     if id=1 then
                     begin
                          writeln(otl,'iro0=',iro0:10:5);
                          writeln(otl,'r=',r:10:5);
                          writeln(otl,'teta0=',iteta0:10:5);
                          writeln(otl,'rot=',rot:10:5);
                          writeln(otl,'tetat=',tetat:10:5);
                          writeln(otl,'massa=',massa:10:5);
                          writeln(otl,'vr=',vr:10:5);
                    end;
                     if ev>1 then exit;
                     if bro0 then
                     begin
                          dro:=0.000001;
                          mpdq.tpvd(iro0+dro,iteta0,rot,tetat,r,vr,massa,et[id],pia,rog1,tetag1,t0,ev);
                          xp:=rog1*sin(tetag1*pi/180);
                          yp:=rog1*cos(tetag1*pi/180);

                          a^[i*2-1,jb+1]:=(xp-x0)/dro;
                          a^[i*2,jb+1]:=(yp-y0)/(dro);
                          jb:=jb+1;
                     end;

                     if bteta0 then
                     begin
                          dro:=0.000001*pi/180;
                          mpdq.tpvd(iro0,iteta0+dro,rot,tetat,r,vr,massa,et[id],pia,rog1,tetag1,t0,ev);
                          xp:=rog1*sin(tetag1*pi/180);
                          yp:=rog1*cos(tetag1*pi/180);

                          a^[i*2-1,jb+1]:=(xp-x0)/dro;
                          a^[i*2,jb+1]:=(yp-y0)/(dro);
                          jb:=jb+1;
                     end;
                     if brot then
                     begin
                          dro:=0.000001;
                          mpdq.tpvd(iro0,iteta0,rot+dro,tetat,r,vr,massa,et[id],pia,rog1,tetag1,t0,ev);
                          xp:=rog1*sin(tetag1*pi/180);
                          yp:=rog1*cos(tetag1*pi/180);

                          a^[i*2-1,jb+1]:=(xp-x0)/dro;
                          a^[i*2,jb+1]:=(yp-y0)/(dro);
                          jb:=jb+1;
                     end;
                     if btetat then
                     begin
                          dro:=0.000001;
                          mpdq.tpvd(iro0,iteta0,rot,tetat+dro,r,vr,massa,et[id],pia,rog1,tetag1,t0,ev);
                          xp:=rog1*sin(tetag1*pi/180);
                          yp:=rog1*cos(tetag1*pi/180);

                          a^[i*2-1,jb+1]:=(xp-x0)/dro;
                          a^[i*2,jb+1]:=(yp-y0)/(dro);
                          jb:=jb+1;
                     end;
                     if br then
                     begin
                          dro:=0.000001;
                          mpdq.tpvd(iro0,iteta0,rot,tetat,r+dro,vr,massa,et[id],pia,rog1,tetag1,t0,ev);
                          {mpdq.tpvd(iro0,iteta0,rot,tetat,r-dro/2,vr,massa,et[id],pia,rog2,tetag1,t0,ev);}

                          xp:=rog1*sin(tetag1*pi/180);
                          yp:=rog1*cos(tetag1*pi/180);

                          a^[i*2-1,jb+1]:=(xp-x0)/dro;
                          a^[i*2,jb+1]:=(yp-y0)/(dro);
                          {writeln(rsx,et[i],' ',a^[i*2-1,jb+1],' ',a^[i*2,jb+1]);}
                          jb:=jb+1;
                     end;


                     if bvr then
                     begin
                          dro:=0.000001;
                          mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr+dro,massa,et[id],pia,rog1,tetag1,t0,ev);
                          xp:=rog1*sin(tetag1*pi/180);
                          yp:=rog1*cos(tetag1*pi/180);

                          a^[i*2-1,jb+1]:=(xp-x0)/dro;
                          a^[i*2,jb+1]:=(yp-y0)/(dro);
                          jb:=jb+1;
                     end;
                     if bmassa then
                     begin
                          dro:=0.000001;
                          mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa+dro,et[id],pia,rog1,tetag1,t0,ev);
                          xp:=rog1*sin(tetag1*pi/180);
                          yp:=rog1*cos(tetag1*pi/180);

                          a^[i*2-1,jb+1]:=(xp-x0)/dro;
                          a^[i*2,jb+1]:=(yp-y0)/(dro);
                          jb:=jb+1;
                     end;

                     if bpia then
                     begin
                          dro:=0.000001;
                          mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et[id],pia+dro,rog1,tetag1,t0,ev);
                          xp:=rog1*sin(tetag1*pi/180);
                          yp:=rog1*cos(tetag1*pi/180);

                          a^[i*2-1,jb+1]:=(xp-x0)/dro;
                          a^[i*2,jb+1]:=(yp-y0)/(dro);
                          jb:=jb+1;
                     end;
                     a^[i*2-1,jb+1]:=x0-x[i];
                     a^[i*2,jb+1]:=y0-y[i];
                     {if (abs(a^[i*2-1,jb+1])>0.5) or (abs(a^[i*2,jb+1])>0.5) then i:=i-1;}
                     {writeln(dt0,et[i]:10:5,' ', rog:10:5,' ',tetag:10:5)}
                end;
                {halt(0);}
                {close(rsx);}
                ptlst(a,xf,cv,i*2,jb);
                sgmnk(a,cv,xf,sgm,jb,i*2);
                     jb:=0;
                     if bro0 then
                     begin
                          jb:=jb+1;
                          xf0[jb]:=iro0;
                     end;

                     if bteta0 then
                     begin
                          jb:=jb+1;
                          xf0[jb]:=iteta0;
                     end;
                     if brot then
                     begin
                          jb:=jb+1;
                          xf0[jb]:=rot;
                     end;
                     if btetat then
                     begin
                          jb:=jb+1;
                          xf0[jb]:=tetat;
                     end;
                     if br then
                     begin
                          jb:=jb+1;
                          xf0[jb]:=r;
                     end;

                     if bvr then
                     begin
                          jb:=jb+1;;
                          xf0[jb]:=vr;
                     end;

                     if bmassa then
                     begin
                          jb:=jb+1;
                          xf0[jb]:=massa;
                     end;

                     if bpia then
                     begin
                          jb:=jb+1;
                          xf0[jb]:=pia;
                     end;


                vyh:=true;
                kvyh:=false;
                for i:=1 to jb do
                begin
                     vyh:=vyh and (abs((xf0[i]-xf2[i])/xf0[i])<1e-5);
                     kvyh:=kvyh or (abs((xf0[i]-xf1[i])/xf0[i])>1e-8);
                end;

                if vyh and (kvyh) and (lp>0.0001) then
                begin
                     {lp:=lp/2;}
                     for i:=1 to 7 do xf1[i]:=0;
                     for i:=1 to 7 do xf2[i]:=0;
                     for i:=1 to 7 do xf0[i]:=0;
                end;
                xf2:=xf1;
                xf1:=xf0;

                     jb:=0;
                     vyh:=true;

                     if bro0 then
                     begin
                          jb:=jb+1;
                          iro0:=iro0-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<0.001);
                     end;

                     if bteta0 then
                     begin
                          jb:=jb+1;
                          iteta0:=iteta0-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<0.001);
                     end;
                     if brot then
                     begin

                          jb:=jb+1;
                          rot:=rot-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<1e-6);
                     end;
                     if btetat then
                     begin
                          jb:=jb+1;
                          tetat:=tetat-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb]/tetat)<0.001);

                     end;
                     if br then
                     begin
                          jb:=jb+1;
                          r:=r-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb]/r)<0.001);

                     end;

                     if bvr then
                     begin
                          jb:=jb+1;;
                          vr:=vr-lp*xf[jb];

                          vyh:=vyh and (abs(xf[jb]/vr)<0.001);

                     end;

                     if bmassa then
                     begin
                          jb:=jb+1;
                          massa:=massa-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb]/massa)<0.001);

                     end;

                     if bpia then
                     begin
                          jb:=jb+1;
                          pia:=pia-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb]/pia)<0.001);
                     end;
                fd:=fd+1;

                if fd mod 1 =0 then
                begin
                gotoxy(1,1);
                writeln;
                if bro0 then pst:='     var'
                        else pst:=' not var'  ;
                writeln('   ro0=',iro0:10:5,pst);
                if bteta0 then pst:='     var'
                        else pst:=' not var'  ;
                writeln(' teta0=',iteta0:10:5,pst);
                if brot then pst:='     var'
                        else pst:=' not var'  ;
                writeln('   rot=',rot:20:10,pst);
                if btetat then pst:='     var'
                        else pst:=' not var'  ;
                writeln(' tetat=',tetat:10:5,pst);
                if br then pst:='     var'
                        else pst:=' not var'  ;
                writeln('krivizn=',r:10:5,pst);
                if bvr then pst:='     var'
                        else pst:=' not var'  ;
                writeln('     vr=',vr:10:5,pst);
                if bmassa then pst:='     var'
                         else pst:=' not var'  ;
                writeln('  massa=',massa:10:5,pst);
                if bpia then pst:='     var'
                         else pst:=' not var'  ;
                writeln('     pi=',pia:10:5,pst);


                writeln('iteracija=',fd,' kol-vo urav=',i,' iz=',iz,' lp=',lp:13:10 );
                if biskl then wvh;
                {if fd mod 10 =0 then
                begin
                     ch:=readkey;
                     if ch='q' then halt(0);
                end;}
                {vydres(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt);;}
                if biskl then vydres(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt,sgm,bro0,bteta0,brot,btetat,br
                ,bvr,bmassa,bpia
                ,nro0,sro0,nteta0,steta0,nrot,srot,ntetat,stetat,oevr,oevt);;

                end;
          until vyh;
          {vydres(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt);;;;}

          {nevotap(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt);;;;}
          if biskl then maxnev(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt,oevr,oevt);
          if biskl then vydres(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt,sgm,bro0,bteta0,brot,btetat,br,bvr
          ,bmassa,bpia,nro0,sro0,nteta0,steta0,nrot,srot,ntetat,stetat,oevr,oevt);;

          if biskl then iskl
                   else ew:=false;
          if biskl then wvh;

          until not(ew);
          if biskl then vydres2(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt,oevr,oevt);
          end;

end;*)
(*procedure vydryd;
var i,j,km,nr,mi,i2,nt2,nt3:integer;
et2,teta2,ro2,x2,y2:dst;
ref:string[6];
t0,min,max,nsr1,nsr2:extended;
rs,rs2:text;
begin
     new(arot);new(atetat);new(sarot);new(satetat);
     new(x3);new(y3);new(et3);new(ro3);new(teta3);
     new(aref3);
     nrd:=0;
     nt2:=0;
     assign(rs,'ryd.res');
     rewrite(rs);
     assign(rs2,'ryd2.res');
     rewrite(rs2);
     for j:=1 to nt do writeln(rs2,j:4,' ',aref^[j]:6,' ',et[j]:8:3,' ',ro[j]:10:5);
     flush(rs2);

     for i:=1 to nt do

     begin
          uge[i]:=false;
     end;
     nr:=0;
     for i:=1 to nt do
     begin
          km:=1;
          ref:=aref^[i];
          mi:=-1;
          for j:=1 to nr do
          begin
               if ref=aref2^[j] then mi:=j;
          end;
          if (mi<0) and (length(ref)>=6) then
          begin
               nr:=nr+1;
               aref2^[nr]:=ref;
               km:=1;
               et2[km]:=et[i];
               teta2[km]:=teta[i];
               ro2[km]:=ro[i];
               x2[km]:=x[i];
               y2[km]:=y[i];
               t0:=et2[km];
               min:=t0;
               max:=t0;

               for j:=i+1 to nt do
               begin
                    if ref=aref^[j] then
                    begin
                         km:=km+1;
                         et2[km]:=et[J];
                         teta2[km]:=teta[j];
                         ro2[km]:=ro[j];
                         x2[km]:=x[j];
                         y2[km]:=y[j];
                         t0:=t0+et2[km];
                         if min>et2[km] then min:=et2[km];
                         if max<et2[km] then max:=et2[km];
                    end;
               end;
               if (km>=4) and (max-min>3.5) then
               begin
                    t0:=t0/km;
                    for j:=1 to km do
                    begin
                         a^[j,1]:=et2[j]-t0;
                         a^[j,2]:=1;
                         a^[j,3]:=ro2[j];
                         at^[j,1]:=et2[j]-t0;
                         at^[j,2]:=1;
                         at^[j,3]:=teta2[j];
                    end;
                    ptlst(a,xf,cv,km,2);
                    sgmnk(a,cv,xf,sgm,2,km);

                    ptlst(at,xft,cv,km,2);
                    sgmnk(at,cv,xft,sgmt,2,km);

                    writeln(rs,ref,' ',km:3,' ',t0:8:3,' ',max-min:8:3,' ',sgm[3]:10:5,' ',sgmt[3]:10:5,
                    ' ',xf[2]:10:5,' ',xft[2]:10:5,
                    ' ',xf[1]:10:5,' ',xft[1]:10:5,
                    ' ',sgm[2]:10:5,' ',sgmt[2]:10:5,
                    ' ',sgm[1]:10:5,' ',sgmt[1]:10:5,
                    ' ',xf[1]/sgm[1]:5:1,' ',xft[1]/sgmt[1]:5:1);
                    nsr1:=xf[1]/sgm[1];
                    nsr2:=xft[1]/sgmt[1];
                    if nsr2>nsr1 then nsr1:=nsr2;
                    if nsr1>3 then
                    begin
                          nrd:=nrd+1;
                          nt2:=nt2+1;
                          arot^[nrd]:=xf[1];
                          atetat^[nrd]:=xft[1];
                          sarot^[nrd]:=sgm[1];
                          satetat^[nrd]:=sgmt[1];

                          et3^[nt2]:=t0;
                          ro3^[nt2]:=xf[2];
                          teta3^[nt2]:=xft[2];
                          x3^[nt2]:=ro[nt2]*sin(teta[nt2]*pi/180);
                          y3^[nt2]:=ro[nt2]*cos(teta[nt2]*pi/180);
                          aref3^[nt2]:=ref;
                    end;


               end;
          end;
     end;
                          writeln(rs2);
                          for j:=1 to nt2 do writeln(rs2,j:4,' ',aref3^[j]:6,' ',et3^[j]:8:3,' ',ro3^[j]:10:5);
                          flush(rs2);
     nt3:=nt2;
     mi:=-1;
     for i:=1 to nt do
     begin
          mi:=-1;
          for j:=1 to nt3 do
          begin
               if aref^[i]=aref3^[j] then mi:=j;
          end;
          if mi<0 then
          begin
               nt2:=nt2+1;
               et3^[nt2]:=et[i];
               ro3^[nt2]:=ro[i];
               teta3^[nt2]:=teta[i];
               x3^[nt2]:=x[i];
               y3^[nt2]:=y[i];
               aref3^[nt2]:=aref^[i];

          end;
     end;
                          writeln(rs2);
                          for j:=1 to nt2 do writeln(rs2,j:4,' ',aref3^[j]:6,' ',et3^[j]:8:3,' ',ro3^[j]:10:5);
                          flush(rs2);
    for j:=1 to nt2 do
    begin
         et[j]:=et3^[j];
         ro[j]:=ro3^[j];
         teta[j]:=teta3^[j];
         x[j]:=x3^[j];
         y[j]:=y3^[j];
         aref^[j]:=aref3^[j];
    end;
    nt:=nt2;


     close(rs);
     close(rs2);
end;*)

procedure ishdata(pst:string);
var pst2,pst3,pst4:string;
p,ln,is,er,nhip:integer;
cha,chz,ch0,ch9:char;

srh,srm:string;
tcat5,tcat2:string[5];
tm1,tm2:extended;
begin
     cha:='A';
     chz:='z';
     ch0:='0';
     ch9:='9';
     assign(dt,PST);
     reset(dt);
     nt:=0;
     sro:=0;
     read(dt,pst12);
     swds:=copy(pst12,2,7)+copy(pst12,10,3);
     if swds[6]='S' then swds[6]:='-' else swds[6]:='+';

     srh:=copy(pst12,2,2);
     val(srh,rh,er);
     erval(srh,er);

     srh:=copy(pst12,4,3);
     val(srh,rm,er);
     erval(srh,er);
     rm:=rm/10.0;

     srh:=pst12[8]+pst12[10];
     val(srh,dg,er);
     erval(srh,er);

     srm:=copy(pst12,11,2);
     val(srh,ds,er);
     erval(srh,er);

     if pst12[7]='S' then
     begin
          dg:=-dg;
          ds:=-ds;
     end;

     ra:=(rh+rm/60)*15*pi/180;
     dc:=(dg+ds/60)*pi/180;

     ra:=sin(ra);
     dc:=1/cos(dc);
     read(dt,ch,comp10);
     {writeln('comp10=',comp10);
     readln;}
     read(dt,sorb1,sorb2,sorb3,sorb4,sorb5,m1,m2);
     bprd:=true;
     if btorb then
     begin
          bprd:=false;

          if sorb1>10 then
          begin
                if sorb2/sorb1<0.3 then
                begin
                      bprd:=true;
                end;
          end;
     end;
     if not(bprd) then
     begin
          close(dt);
          exit;
     end;

     read(dt,nhip);
     if nhip>0 then
     begin
          read(dt,pia);
          pia:=pia/1000;
     end else
     begin
          pia:=1;
     end;
     readln(dt);
     r26:=false;
     while not(eof(dt)) do
     begin
          read(dt,tcat5,ch,tcomp10);
          tcat2:=copy(tcat5,3,2);
          r26:=r26 or (tcat2='26');
                nt:=nt+1;
                read(dt,et^[nt],teta^[nt],ro^[nt],tm1,tm2,pst2);
                aref^[nt]:=copy(pst2,2,8);
                if tcat5='cd26b' then
                begin
                     aref^[nt]:='ccd26ref';
                end;
                if tcat5='ph26c' then
                begin
                     aref^[nt]:='pht26ref';
                end;
                pia:=1;


                {fed:=copy(pst,1,8);

                fed:=otnc(fed);
                val(fed,et[nt],ercode);
                erval(fed,ercode);

                fed:=copy(pst,13,7);
                fed:=otnc(fed);
                val(fed,teta[nt],ercode);
                erval(fed,ercode);

                fed:=copy(pst,22,6);
                fed:=otnc(fed);
                val(fed,ro[nt],ercode);
                erval(fed,ercode);}


                {if nt<=t3 then}

                begin
                     teta^[nt]:=teta^[nt]+0.00556*ra*(2000-et^[nt])*dc;
                end;

                sro:=sro+ro^[nt];
                x^[nt]:=ro^[nt]*sin(teta^[nt]*pi/180);
                y^[nt]:=ro^[nt]*cos(teta^[nt]*pi/180);

                {readln(dt,pst2);
                p:=pos(':',pst2);
                ln:=length(pst2);
                pst3:='';
                is:=0;
                while (is<3) and (p<ln) do
                begin
                     p:=p+1;
                     if (cha<=pst2[p]) and (pst2[p]<=chz) then
                     begin
                          is:=is+1;
                          pst3:=pst3+pst2[p];
                     end;
                end;
                p:=ln;
                pst4:='';
                is:=0;
                while (is<3) and (p>1) do
                begin
                     if (ch0<=pst2[p]) and (pst2[p]<=ch9) then
                     begin
                          is:=is+1;
                          pst4:=pst2[p]+pst4;
                     end;
                     p:=p-1;

                end;
                pst3:=pst3+pst4;
                aref^[nt]:=pst3;}
                p:=0;
                readln(dt);
     end;
     t4:=1;

     for i:=1 to nt do et^[i]:=et^[i]-t0;
     nnt:=nt;
     t1:=t4;
     close(dt);
     writeln('nt=',nt);
     if (pst12='S17066N0.039') and (pos('Ba,Bb',comp10)>0) then
     begin
          teta^[1]:=teta^[1]-180;
          x^[1]:=ro^[1]*sin(teta^[1]*pi/180);
          y^[1]:=ro^[1]*cos(teta^[1]*pi/180);

     end;

     {if nt>500 then nt:=1;}

end;

procedure prf;
begin
      assign(dt,'vyorbit.prf');
      reset(dt);

      readln(dt,ch);
      if ch='T' then bro0:=true;
      if ch='F' then bro0:=false;
      if (ch<>'T') and (ch<>'F') then
      begin
            writeln(' Неверный формат файла vyorbit.prf. Cтрока 1.');
            halt(0);
      end;

      readln(dt,ch);
      if ch='T' then bteta0:=true;
      if ch='F' then bteta0:=false;
      if (ch<>'T') and (ch<>'F') then
      begin
            writeln(' Неверный формат файла vyorbit.prf. Cтрока 2.');
            halt(0);
      end;

      readln(dt,ch);
      if ch='T' then brot:=true;
      if ch='F' then brot:=false;
      if (ch<>'T') and (ch<>'F') then
      begin
            writeln(' Неверный формат файла vyorbit.prf. Cтрока 3.');
            halt(0);
      end;

      readln(dt,ch);
      if ch='T' then btetat:=true;
      if ch='F' then btetat:=false;
      if (ch<>'T') and (ch<>'F') then
      begin
            writeln(' Неверный формат файла vyorbit.prf. Cтрока 4.');
            halt(0);
      end;

      readln(dt,ch);
      if ch='T' then br:=true;
      if ch='F' then br:=false;
      if (ch<>'T') and (ch<>'F') then
      begin
            writeln(' Неверный формат файла vyorbit.prf. Cтрока 5.');
            halt(0);
      end;

      readln(dt,ch);
      if ch='T' then bvr:=true;
      if ch='F' then bvr:=false;
      if (ch<>'T') and (ch<>'F') then
      begin
            writeln(' Неверный формат файла vyorbit.prf. Cтрока 6.');
            halt(0);
      end;

      readln(dt,ch);
      if ch='T' then bmassa:=true;
      if ch='F' then bmassa:=false;
      if (ch<>'T') and (ch<>'F') then
      begin
            writeln(' Неверный формат файла vyorbit.prf. Cтрока 7.');
            halt(0);
      end;

      readln(dt,ch);
      if ch='T' then bpia:=true;
      if ch='F' then bpia:=false;
      if (ch<>'T') and (ch<>'F') then
      begin
            writeln(' Неверный формат файла vyorbit.prf. Cтрока 8.');
            halt(0);
      end;

      readln(dt,l);
      le:=l;
      readln(dt,rh,rm,dg,ds);
      readln(dt,t0);
      readln(dt,r);
      readln(dt,vr);
      readln(dt,massa);
      readln(dt,pia);

      readln(dt,ch);
      if ch='T' then bnkri:=true;
      if ch='F' then bnkri:=false;
      if (ch<>'T') and (ch<>'F') then
      begin
            writeln(' Неверный формат файла vyorbit.prf. Cтрока 16.');
            halt(0);
      end;

      readln(dt,ch);
      if ch='T' then bnvr:=true;
      if ch='F' then bnvr:=false;
      if (ch<>'T') and (ch<>'F') then
      begin
            writeln(' Неверный формат файла vyorbit.prf. Cтрока 17.');
            halt(0);
      end;

      readln(dt,ch);
      if ch='T' then bprr:=true;
      if ch='F' then bprr:=false;
      if (ch<>'T') and (ch<>'F') then
      begin
            writeln(' Неверный формат файла vyorbit.prf. Cтрока 18.');
            halt(0);
      end;

      readln(dt,ch);
      if ch='T' then bos:=true;
      if ch='F' then bos:=false;
      if (ch<>'T') and (ch<>'F') then
      begin
            writeln(' Неверный формат файла vyorbit.prf. Cтрока 19.');
            halt(0);
      end;
      close(dt);
end;

(*procedure genish;

var i,j:integer;
    dt:text;
    rop,tetap,xp,yp,r:extended;
begin

     assign(dt,'genish.dt');
     rewrite(dt);
     writeln(dt,t100);

     for j:=1 to t100 do
     begin
          massa:=massa0+smassa*(random+random+random+
                                random+random+random+
                                random+random+random+
                                random+random+random-6);
          pia:=pi0+spia*(random+random+random+
                         random+random+random+
                         random+random+random+
                         random+random+random-6);
          r:=kriv0+skri*(random+random+random+
                         random+random+random+
                         random+random+random+
                         random+random+random-6);
          write(dt,massa:10:4,' ',pia:12:7,' ',nt,' ');
          for i:=1 to nt do
          begin
               rop:=xro0^[i]+oevr0*(random+random+random+
                                      random+random+random+
                                      random+random+random+
                                      random+random+random-6);

               tetap:=xteta0^[i]+oevt0*(random+random+random+
                                          random+random+random+
                                          random+random+random+
                                          random+random+random-6);
               write(dt,' ',et[i]:10:4,' ',rop:8:4,' ',tetap:8:4);

          end;
          writeln(dt);
     end;
     close(dt);
end;

procedure genish2;
begin
          read(dtp,massa,pia,nt);
          for i:=1 to nt do
          begin
               read(dtp,et[i],ro[i],teta[i]);
               x[i]:=ro[i]*sin(teta[i]*pi/180);
               y[i]:=ro[i]*cos(teta[i]*pi/180);
          end;
          readln(dtp);
end;*)

procedure rvhd;
var dt:text;
begin
      assign(dt,'rvhd.dt');
      reset(dt);
      readln(dt,iz,iv);

      readln(dt,abp0);
      readln(dt,ev0);
      readln(dt,ai0);
      readln(dt,omm0);
      readln(dt,omb0);
      readln(dt,th0);
      readln(dt,prd0);

      readln(dt,sabp0);
      readln(dt,sev0);
      readln(dt,sai0);
      readln(dt,somm0);
      readln(dt,somb0);
      readln(dt,sth0);
      readln(dt,sprd0);
      for i:=2000 to t2050 do
      begin
           readln(dt,roe^[i],tetae^[i],sroe^[i],statae^[i]);
      end;

      close(dt);
      assign(dtp,'genish.dt');
      reset(dtp);
      readln(dtp);
      for i:=1 to iz do readln(dtp);
end;
procedure wvhd;
var dt:text;
begin
      assign(dt,'rvhd.dt');
      rewrite(dt);
      writeln(dt,iz,' ',iv);

      writeln(dt,abp0);
      writeln(dt,ev0);
      writeln(dt,ai0);
      writeln(dt,omm0);
      writeln(dt,omb0);
      writeln(dt,th0);
      writeln(dt,prd0);

      writeln(dt,sabp0);
      writeln(dt,sev0);
      writeln(dt,sai0);
      writeln(dt,somm0);
      writeln(dt,somb0);
      writeln(dt,sth0);
      writeln(dt,sprd0);
      for i:=2000 to t2050 do
      begin
           writeln(dt,roe^[i],' ',tetae^[i],' ',sroe^[i],' ',statae^[i]);
      end;
      close(dt);
end;

procedure wvhd0;
var dt:text;
begin
      assign(dt,'rvhd.dt');
      rewrite(dt);
      writeln(dt,'0 0');

      writeln(dt,abp0);
      writeln(dt,ev0);
      writeln(dt,ai0);
      writeln(dt,omm0);
      writeln(dt,omb0);
      writeln(dt,th0);
      writeln(dt,prd0);

      writeln(dt,'0');
      writeln(dt,'0');
      writeln(dt,'0');
      writeln(dt,'0');
      writeln(dt,'0');
      writeln(dt,'0');
      writeln(dt,'0');
      for i:=2000 to t2050 do
      begin
           writeln(dt,roe^[i],' ',tetae^[i],' 0 0');
      end;
      close(dt);
end;

(*procedure podgotovka(gh:boolean);
begin
     if gh then
     begin
           gen0(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt,oevr,oevt);
           genish;
           wvhd0;
     end;
end;*)

function dl2(iro0,iteta0,rot,tetat,r,vr,massa,pia:extended):extended;
var
ti0:das;
dro,dte,raz,et2,rog,tetag,sum,h,x,y,xp,yp,dr:extended;
id,ig,km,dk,mig,nt:integer;
ep,ro,teta,max,ete,etb,pte,max2,sum2:extended;
dt:text;
begin
     km:=0;
     for ig:=1 to nnt do
     begin
          if gbop^[ig]=0 then
          begin
               km:=km+1;
               aex^[km]:=et^[ig];
               aex2^[km]:=0;
               aex3^[km]:=0;
          end;
     end;


      mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,vr,massa,0,pia,ti0);
      et2:=0;
      nt:=1000;
      sum:=0;
      ete:=aex^[km];
      etb:=aex^[1];

      h:=(ete-etb)/nt;
      sum2:=0;
      for id:=0 to nt do
      begin
           et2:=etb+id*h;
           mpdq.rtti(ti0,et2,rog,tetag);
           x:=rog*cos(tetag*pi/180);
           y:=rog*sin(tetag*pi/180);
           if id<>0 then
           begin
                dr:=sqrt(sqr(x-xp)+sqr(y-yp));
                sum:=sum+dr;
                dte:=abs(tetag-pte);
                ogr2(dte);
                sum2:=sum2+dte;
                for ig:=1 to km-1 do
                begin
                     if (aex^[ig]<=et2) and (et2<=aex^[ig+1]) then
                     begin
                          aex2^[ig]:=aex2^[ig]+dr;
                          aex3^[ig]:=aex3^[ig]+dte;
                     end;
                end;
           end;
           xp:=x;
           yp:=y;
           pte:=tetag;
      end;
      max:=-1e39;
      max2:=-1e39;
      for ig:=1 to km-1 do
      begin
           if aex2^[ig]>max then
           begin
                max:=aex2^[ig];
                mig:=ig;
           end;
           if aex3^[ig]>max2 then
           begin
                max2:=aex3^[ig];
           end;
      end;
      {dl2:=sum-max;}
      sum:=sum-max;
      sum2:=sum2-max2;
      if sum<0.000001 then sum:=0.000001;
      {dl2:=max/sum;}
      if sum2<0.000001 then sum2:=0.000001;

      dl2:=max2/sum2;


end;
procedure telcir3(var ib:integer;
                  var siro0,siteta0,srot,stetat,sr,sgf0,sgprd,sgnak
                      ,gonak,gonak2,sgonak,sgonak2,pdt,pdx,pdy:extended;
                  var ugep,bbiskl:boolean);
var ix,is,ig:integer;
tdl2,lf0,te1,te2,raz1,raz2:extended;
bgmi,bp180:boolean;
ti0:das;
begin
               siro0:=iro0;
               siteta0:=iteta0;
               srot:=rot;
               stetat:=tetat;
               sr:=r;
               for ix:=1 to 6 do svsgm[ix]:=sgm[ix];
               if sgf0<10 then  le:=0.01
                          else  le:=1e-6;

               for is:=1 to nnt do gbop2^[is]:=gbop^[is];
               ciriteracii(bbiskl);
               tdl2:=dl2(iro0,iteta0,rot,tetat,r,vr,massa,pia);

               if tdl2>tt10 then gf0:=1e39;
               if gf0>5e38 then for is:=1 to nnt do gbop^[is]:=gbop2^[is];

               gonak:=abs(90-gnak*180/pi);
               gonak2:=abs(180-gnak*180/pi);
               if gonak2<gonak then gonak:=gonak2;
               if gonak<0.01 then gonak:=0.01;
               if sgprd<0.01 then sgprd:=0.01;
               if sgf0<0.0001 then sgf0:=0.0001;

               if ((gmi=ib) and ugep)
                  or (gbop^[ib]=2)
                  or ((gf0/sgf0)>2)
                  or (gf0>5e38)
                  or ((gprd/sgprd)>2.5)
                  or ((sgonak/gonak)>4)
                  then
               begin
                    iro0:=siro0;
                    iteta0:=siteta0;
                    rot:=srot;
                    tetat:=stetat;
                    r:=sr;
                    for ix:=1 to 6 do sgm[ix]:=svsgm[ix];
                    bp180:=false;
                    if abs(m2-m1)<1 then
                    begin

                          if not((r<0) or (r>3000)) and not((massa<0) or (massa>3000)) then
                          begin
                                mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,vr,massa,0,pia,ti0);
                                mpdq.rtti(ti0,et^[ib],rog,tetag);
                                te1:=teta^[ib];
                                te2:=te1-180;
                                if te2<0 then te2:=te2+360;
                                raz1:=(tetag-te1);
                                ogr2(raz1);
                                RAZ1:=abs(raz1*pi*ro^[ib]/180);

                                raz2:=(tetag-te2);
                                ogr2(raz2);
                                RAZ2:=abs(raz2*pi*ro^[ib]/180);

                                bp180:=raz1>10*raz2;

                          end;
                    end;
                    if (abs(m2-m1)<1) and bp180 then
                    begin
                         pdt:=teta^[ib];
                         pdx:=x^[ib];
                         pdy:=y^[ib];
                         teta^[ib]:=teta^[ib]-180;
                         if teta^[ib]<0 then teta^[ib]:=teta^[ib]+360;
                         x^[ib]:=ro^[ib]*sin(teta^[ib]*pi/180);
                         y^[ib]:=ro^[ib]*cos(teta^[ib]*pi/180);
                         gbop^[ib]:=0;
                        if sgf0<10 then  le:=0.01
                                   else  le:=1e-6;

                         for is:=1 to nnt do gbop2^[is]:=gbop^[is];
                         ciriteracii(bbiskl);
                         tdl2:=dl2(iro0,iteta0,rot,tetat,r,vr,massa,pia);
                         if tdl2>tt10 then gf0:=1e39;
                         if gf0>5e38 then for is:=1 to nnt do gbop^[is]:=gbop2^[is];


                         gonak:=abs(90-gnak*180/pi);
                         gonak2:=abs(180-gnak*180/pi);
                         if gonak2<gonak then gonak:=gonak2;
                         if gonak<0.01 then gonak:=0.01;
                         if sgprd<0.01 then sgprd:=0.01;
                         if sgf0<0.0001 then sgf0:=0.0001;

                        if
                           ((gmi=ib) and ugep)
                           or (gbop^[ib]=2)
                           or ((gf0/sgf0)>2)
                           or ((gprd/sgprd)>2.5)
                           or ((sgonak/gonak)>4)
                           or (gf0>5e38)

                        then
                        begin
                             teta^[ib]:=pdt;
                             x^[ib]:=pdx;
                             y^[ib]:=pdy;
                             iro0:=siro0;
                             iteta0:=siteta0;
                             rot:=srot;
                             tetat:=stetat;
                             r:=sr;
                             for ix:=1 to 6 do sgm[ix]:=svsgm[ix];

                             if gbop^[ib]=2 then gbop^[ib]:=2;
                             if (gmi=ib) then gbop^[ib]:=3;
                             if ((gf0/sgf0)>2) then gbop^[ib]:=4;
                             if ((gprd/sgprd)>2.5) then gbop^[ib]:=5;
                             if ((sgonak/gonak)>4) then gbop^[ib]:=6;
                             if gbop^[ib]=0 then  gbop^[ib]:=1;


                        end else
                        begin
                             if sgf0>gf0 then sgf0:=gf0;
                             if sgprd>gprd then sgprd:=gprd;
                             if sgonak<gonak then sgonak:=gonak;
                             ugep:=true;
                             if sgf0<0.05 then sgf0:=0.05;
                        end;
                    end else
                    begin
                             if gbop^[ib]=2 then gbop^[ib]:=2;
                             if (gmi=ib) then gbop^[ib]:=3;
                             if ((gf0/sgf0)>2) then gbop^[ib]:=4;
                             if ((gprd/sgprd)>2.5) then gbop^[ib]:=5;
                             if ((sgonak/gonak)>4) then gbop^[ib]:=6;
                             if gbop^[ib]=0 then  gbop^[ib]:=1;

                    end;
               end else
               begin
                    if sgf0>gf0 then sgf0:=gf0;
                    if sgprd>gprd then sgprd:=gprd;
                    if sgonak<gonak then sgonak:=gonak;
                    if sgf0<0.05 then sgf0:=0.05;

                    ugep:=true;

               end;


end;
procedure telcir2(var ib:integer;
                  var siro0,siteta0,srot,stetat,sr,sgf0,sgprd,sgnak
                      ,gonak,gonak2,sgonak,sgonak2,pdt,pdx,pdy:extended;
                  var ugep,bbiskl:boolean);
var ix,is,ig:integer;
tdl2,lf0:extended;
bgmi:boolean;
begin
          if gbop3^[ib]=1 then
          begin

               gbop^[ib]:=0;
               lf0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
               bzc:=false;
               ngmi:=ngmi+1;
               agmi^[ngmi]:=nevgmi;
               bgmi:=false;
               for ig:=1 to ngmi do bgmi:=bgmi or (agmi^[ig]=ib);
               if sgf0<0.0001 then sgf0:=0.0001;

               if bgmi or ((lf0/gf0)>1.2) then
               begin
                       {telcir3}
                    ngmi:=0;
                  telcir3(ib,
                          siro0,siteta0,srot,stetat,sr,sgf0,sgprd,sgnak
                          ,gonak,gonak2,sgonak,sgonak2,pdt,pdx,pdy,
                           ugep,bbiskl);
               end;
          end else
          begin
               ugep:=true;
          end;

end;
procedure predv(bbiskl:boolean);
var siro0,siteta0,srot,stetat,sr:extended;
xsv,ysv,rosv,tetasv,etsv:array[1..3] of extended;
gsv:array[1..3] of integer;
ib,t3sv:integer;
begin
      siro0:=iro0;
      siteta0:=iteta0;
      srot:=rot;
      stetat:=tetat;
      sr:=r;
      for ib:=1 to 3 do
      begin
           xsv[ib]:=x^[ib];
           ysv[ib]:=y^[ib];
           rosv[ib]:=ro^[ib];
           tetasv[ib]:=teta^[ib];
           etsv[ib]:=et^[ib];
           gsv[ib]:=gbop^[ib];

           x^[ib]:=ae9[ib*3-2];
           y^[ib]:=ae9[ib*3-1];
           et^[ib]:=ae9[ib*3];
           ro^[ib]:=sqrt(x^[ib]*x^[ib]+y^[ib]*y^[ib]);

           dxw:=x^[ib];
           dyw:=y^[ib];;
           teta^[ib]:=arctan(dxw/dyw);
           if dyw<0 then teta^[ib]:=pi+teta^[ib];
           if teta^[ib]<0 then teta^[ib]:=teta^[ib]+2*pi;
           teta^[ib]:=teta^[ib]*180/pi;
           gbop^[ib]:=0;
      end;
      t3sv:=nnt;
      nnt:=3;
      t3:=3;
      nt:=3;
      ciriteracii(false);
      nnt:=t3sv;
      t3:=t3sv;
      nt:=t3sv;
      for ib:=1 to 3 do
      begin
           x^[ib]:=xsv[ib];
           y^[ib]:=ysv[ib];
           ro^[ib]:=rosv[ib];
           teta^[ib]:=tetasv[ib];
           et^[ib]:=etsv[ib];
           gbop^[ib]:=gsv[ib];
      end;

      if not(uci) then
      begin
            iro0:=siro0;
            iteta0:=siteta0;
            rot:=srot;
            tetat:=stetat;
            r:=sr;
      end;

end;

{procedure getmxt(iro0,iteta0,rot,tetat,r:extended;var max:extended; var mi:integer);
var fln:string;
    dt:text;
    ep,x,y,ro,teta,rog,tetag:extended;
    ib,km:integer;
    ti0:das;
    vr,massa,pia,rteta,dte:extended;
begin
     pia:=1;
     massa:=1;
     vr:=1;
     mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,vr,massa,0,pia,ti0);
     km:=0;
     max:=-1e39;
     for ib:=1 to nnt do
     begin
          if  gbop^[ib]=0 then
          begin
               mpdq.rtti(ti0,et^[ib],rog,tetag);
               km:=km+1;
               if km>1 then
               begin
                    if tetat>0 then
                    begin
                         dte:=tetag-rteta;
                    end else
                    begin
                         dte:=-tetag+rteta;
                    end;
                    km:=km+1;
                    if dte<0 then
                    begin
                         dte:=dte+360;
                    end;
                    if dte>max then
                    begin
                         max:=dte;
                         mi:=ib;
                    end;
               end;
               rteta:=tetag;

          end;
     end;
end;}
procedure ciriteracii3(bbiskl:boolean;
          var ib1,ib2,mi:integer;
          var siro0,siteta0,srot,stetat,sr,sgf0,sgprd,sgnak
              ,gonak,gonak2,sgonak,sgonak2,pdt,pdx,pdy,mxt:extended;
          var ugep:boolean);
var ib:integer;
begin
     ngmi:=0;
     for ib:=1 to nnt do
     begin
          gbop3^[ib]:=gbop^[ib];
     end;
     ib:=0;
     repeat
           ib:=ib+1;
           telcir2(ib,
                   siro0,siteta0,srot,stetat,sr,sgf0,sgprd,sgnak,
                   gonak,gonak2,sgonak,sgonak2,pdt,pdx,pdy,
                   ugep,bbiskl);
     until (ib>=nnt) or (ugep);
     ib1:=ib+1;
     if ib<nnt then
     begin
          ib:=nnt+1;
          ugep:=false;
          while not(ugep) and (ib>=2) do
          begin
               ib:=ib-1;
               telcir2(ib,
                   siro0,siteta0,srot,stetat,sr,sgf0,sgprd,sgnak,
                   gonak,gonak2,sgonak,sgonak2,pdt,pdx,pdy,
                   ugep,bbiskl);
          end;
          ib2:=ib-1;
          for ib:=ib1 to ib2 do
          begin
               if ib=7 then
               begin
                    ib:=7;
               end;
               telcir2(ib,
                   siro0,siteta0,srot,stetat,sr,sgf0,sgprd,sgnak,
                   gonak,gonak2,sgonak,sgonak2,pdt,pdx,pdy,
                   ugep,bbiskl);
          end;

     end;


end;



procedure ciriteracii2(bbiskl:boolean);
var ib,ib1,ib2,mi,ib3,is:integer;
   siro0,siteta0,srot,stetat,sr,sgf0,sgprd,sgnak
   ,gonak,gonak2,sgonak,sgonak2,pdt,pdx,pdy,mxt,tdl2:extended;
   ugep:boolean;

   {diro0,diteta0,drot,dtetat,dr:extended;
   dgf0,dgprd,dgnak,dgonak:extended;}

begin
     t50000:=50000;
     predv(bbiskl);

     for is:=1 to nnt do gbop2^[is]:=gbop^[is];
     ciriteracii(bbiskl);
     tdl2:=dl2(iro0,iteta0,rot,tetat,r,vr,massa,pia);
     if tdl2>tt10 then gf0:=1e39;
     if gf0>5e38 then for is:=1 to nnt do gbop^[is]:=gbop2^[is];


      {diro0:=iro0;
      diteta0:=iteta0;
      drot:=rot;
      dtetat:=tetat;
      dr:=r;}

     t50000:=7500;
     uci2:=uci;
                     if uci2 then
                     begin
                          write(rso,' 8 ');
                     end else
                     begin
                          write(rso,' 0 ');
                    end;

     {if not(posf) and not(uci) then
     begin
          exit;
     end;}
     if gf0<5e38 then
     begin
          sgf0:=gf0;
          sgprd:=gprd;
          sgnak:=gnak;
          sgonak:=abs(90-gnak*180/pi);
          sgonak2:=abs(180-gnak*180/pi);
          if sgonak2<sgonak then sgonak:=sgonak2;
          if sgf0<0.05 then sgf0:=0.05;


     end else
     begin
          sgf0:=gf0;
          sgprd:=1000000;
          sgnak:=pi/4;
          sgonak:=0;

     end;
    ugep:=false;

     {m2:=99.9;}

     {dgf0:=sgf0;
     dgprd:=sgprd;
     dgnak:=sgnak;
     dgonak:=sgonak;}



     ciriteracii3(bbiskl,
                  ib1,ib2,mi,
                  siro0,siteta0,srot,stetat,sr,sgf0,sgprd,sgnak,
                  gonak,gonak2,sgonak,sgonak2,pdt,pdx,pdy,mxt,
                  ugep);

     {getmxt(iro0,iteta0,rot,tetat,r,mxt,mi);
     ib:=0;
     if mxt>82 then
     begin
          for ib3:=1 to mi-1 do
          begin
               aex4^[ib3]:=teta^[ib3];
               teta^[ib3]:=teta^[ib3]-180;
               if teta^[ib3]<0 then teta^[ib3]:=teta^[ib3]+360;
                x^[i]:=ro^[i]*sin(teta^[i]*pi/180);
                y^[i]:=ro^[i]*cos(teta^[i]*pi/180);
          end;
          m2:=99.9;

          iro0:=diro0;
          iteta0:=diteta0;
          rot:=drot;
          tetat:=dtetat;
          r:=dr;
          sgf0:=dgf0;
          sgprd:=dgprd;
          sgnak:=dgnak;
          sgonak:=dgonak;

          for ib:=1 to nnt do
          begin
               gbop^[ib]:=gbop3^[ib];
          end;


          ciriteracii3(bbiskl,
                       ib1,ib2,mi,
                       siro0,siteta0,srot,stetat,sr,sgf0,sgprd,sgnak,
                       gonak,gonak2,sgonak,sgonak2,pdt,pdx,pdy,mxt,
                       ugep);

     end;}
     ib:=0;



end;

procedure ItrIsklB7;
var
ib,is,i2,j2,i3,k3:integer;
dgf0,f0:extended;
begin

                     k3:=0;
                     dgf0:=gf0;

                     repeat
                     i:=1;
                     i2:=-1;
                     while (gbop^[i]>0) and (i<nnt) do
                     begin
                          i:=i+1;
                     end;
                     if (i<nnt) and (i>1) then
                     begin
                          i2:=i-1;
                          while (gbop^[i2]<8) and (i2>1) do
                          begin
                               i2:=i2-1;
                          end;
                          if (i2=1) then
                          begin
                                if (gbop^[1]>=8) then i2:=1 else i2:=-1;
                          end;
                     end else
                     begin
                           if (i=1) and (gbop^[1]>=8) then i2:=1;
                           if (i=nnt) and (gbop^[nnt]>=8) then i2:=nnt;
                     end;
                     j:=nnt;
                     while (gbop^[j]>0) and (j>1) do
                     begin
                          j:=j-1;
                     end;
                     j2:=-1;
                     if (j<nnt) and (j>1) then
                     begin
                          j2:=j+1;
                          while (gbop^[j2]<8) and (j2<nnt) do
                          begin
                               j2:=j2+1;
                          end;
                          if (j2=nnt) then
                          begin
                               if (gbop^[j2]>=8) then j2:=nnt else j2:=-1;
                          end;
                     end else
                     begin
                          if (j=nnt) and (gbop^[j]>=8) then j2:=nnt;
                          if (j=1)   and (gbop^[j]>=8) then j2:=1;
                     end;
                     i3:=-1;
                     if (j2>0) and (i2>0) then
                     begin
                          if abs(et^[j2]-et^[j])<abs(et^[i2]-et^[i]) then
                          begin
                                i3:=j2;
                          end else
                          begin
                               i3:=i2;
                          end;
                     end else
                     begin
                           if j2>0 then i3:=j2;
                           if i2>0 then i3:=i2;
                     end;
                     if i3>0 then
                     begin
                          gbop^[i3]:=0;
                          k3:=k3+1;
                          if l<1e-5 then l:=1e-5;
                          f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
                          if (f0/dgf0)>1.5 then
                          begin
                               riteracii(true);
                               dgf0:=gf0;
                               {writeln('otb');
                               readln;}
                               {halt(0);}
                          end;
                          {ciriteracii(true);}



                     end;


                     until i3<0;

end;

procedure ItrIskl1;
var
ib,is,i2,j2,i3,k3:integer;
dgf0,f0:extended;
begin
                     dgf0:=gf0;
                     k3:=0;
                     repeat

                     i:=1;
                     i2:=-1;
                     while (gbop^[i]>0) and (i<nnt) do
                     begin
                          i:=i+1;
                     end;
                     if (i<nnt) and (i>1) then
                     begin
                          i2:=i-1;
                          while (gbop^[i2]<>1) and (i2>1) do
                          begin
                               i2:=i2-1;
                          end;
                          if (i2=1) then
                          begin
                                if (gbop^[1]=1) then i2:=1 else i2:=-1;
                          end;
                     end else
                     begin
                           if (i=1) and (gbop^[1]=1) then i2:=1;
                           if (i=nnt) and (gbop^[nnt]=1) then i2:=nnt;
                     end;
                     j:=nnt;
                     while (gbop^[j]>0) and (j>1) do
                     begin
                          j:=j-1;
                     end;
                     j2:=-1;
                     if (j<nnt) and (j>1) then
                     begin
                          j2:=j+1;
                          while (gbop^[j2]<>1) and (j2<nnt) do
                          begin
                               j2:=j2+1;
                          end;
                          if (j2=nnt) then
                          begin
                               if (gbop^[j2]=1) then j2:=nnt else j2:=-1;
                          end;
                     end else
                     begin
                          if (j=nnt) and (gbop^[j]=1) then j2:=nnt;
                          if (j=1)   and (gbop^[j]=1) then j2:=1;
                     end;
                     i3:=-1;
                     if (j2>0) and (i2>0) then
                     begin
                          if abs(et^[j2]-et^[j])<abs(et^[i2]-et^[i]) then
                          begin
                                i3:=j2;
                          end else
                          begin
                               i3:=i2;
                          end;
                     end else
                     begin
                           if j2>0 then i3:=j2;
                           if i2>0 then i3:=i2;
                     end;
                     if i3>0 then
                     begin
                          gbop^[i3]:=0;
                          k3:=k3+1;
                          if l<1e-5 then l:=1e-5;
                          f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
                          if (f0/dgf0)>1.1 then
                          begin
                               riteracii(true);
                               dgf0:=gf0;
                               {writeln('otb');
                               readln;}
                               {halt(0);}
                          end;
                          {ciriteracii(true);}



                     end;
                     until i3<0;
                     for i3:=1 to nnt do
                     begin
                          if gbop^[i3]=1 then
                          begin
                                gbop^[i3]:=0;
                                k3:=k3+1;
                                if l<1e-5 then l:=1e-5;
                                f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
                                if (f0/dgf0)>1.1 then
                                begin
                                     riteracii(true);
                                     dgf0:=gf0;
                                end;
                          end;

                     end;

end;
procedure ItrIskl2;
var i3,k3:integer;
dgf0,f0:extended;
begin
                     dgf0:=gf0;
                     k3:=0;
                     for i3:=1 to nnt do
                     begin
                          if gbop^[i3]=2 then
                          begin
                                gbop^[i3]:=0;
                                k3:=k3+1;
                                if l<1e-5 then l:=1e-5;
                                f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
                                if (f0/dgf0)>1.1 then
                                begin
                                     riteracii(true);
                                     dgf0:=gf0;
                                end;
                          end;

                     end;


end;
procedure regular;
var max,min,x1,x2,y1,y2,ras,sras,sras2:extended;
i,j,k,ib,mi:integer;
begin
     max:=-1e39;
     min:=1e39;
     sras:=0;
     mi:=-1;
           for ib:=1 to nnt-1 do
           begin
                RAS:=sqrt(sqr(x^[ib]-x^[ib+1])+sqr(y^[ib]-y^[ib+1]));
                if ras>max then
                begin
                     max:=ras;
                     mi:=ib;
                end;
                if ras<min then
                begin
                     min:=ras;
                end;
                sras:=sras+ras;
           end;
     nnt2:=1;
     sras2:=0;
     k:=0;
     x1:=0;
     y1:=0;
     x2:=0;
     for ib:=1 to nnt do
     begin
          if ib>1 then
          begin
                RAS:=sqrt(sqr(x^[ib]-x^[ib-1])+sqr(y^[ib]-y^[ib-1]));
                sras2:=sras2+ras;
                if sras2>=max then
                begin
                     aex^[nnt2]:=x1/k;
                     aex2^[nnt2]:=y1/k;
                     aex3^[nnt2]:=x2/k;
                     k:=0;
                     nnt2:=nnt2+1;
                     sras2:=0;
                     x1:=0;
                     y1:=0;
                     x2:=0;
                end;
          end;
          x1:=x1+x^[ib];
          y1:=y1+y^[ib];
          x2:=x2+et^[ib];
          k:=k+1;
     end;
     if k>0 then
     begin
          aex^[nnt2]:=x1/k;
          aex2^[nnt2]:=y1/k;
          aex3^[nnt2]:=x2/k;
     end;
     for ib:=1 to nnt2 do
     begin
          x^[ib]:=aex^[ib];
          y^[ib]:=aex2^[ib];
          et^[ib]:=aex3^[ib];
          dxw:=x^[ib];
          dyw:=y^[ib];
          om0:=arctan(dxw/dyw);
          if dyw<0 then om0:=pi+om0;
          if om0<0 then om0:=om0+2*pi;
          teta^[ib]:=om0*180/pi;
          ro^[ib]:=sqrt(dxw*dxw+dyw*dyw);

     end;
     nt:=nnt2;
     nnt:=nnt2;

     k:=0;

end;
procedure vydgbop;
var rs:text;
id:integer;
begin
     assign(rs,'flgbop.res');
     rewrite(rs);
     for id:=1 to nt do
     begin
          writeln(rs,' ',gbop^[id],' ',et^[id]+t0:10:5,' ',ro^[id]:10:5,' ',teta^[id]:10:5);
     end;
     close(rs);
end;
procedure chuzhel3;
var abp,ev,th,prd,aa,bb,ff,gg,comm,somm,comb,somb,ci,si,sum,raz,dte,xa,ya
,xg,yg,sum2,srdro:extended;
ti0:das;
nt2:integer;
ai,omb,omm,ro2014,teta2014:extended;
begin
     prd:=39.4;
     abp:=0.70;
     ev:=0.35;
     ai:=111.6;
     omb:=24.5;
     omm:=354.2 ;
     th:=1999.1-t0;
     getel(swds,comp10,abp,ev,th,prd,ai,omm,omb,ro2014,teta2014);
     if abp>0 then
     begin

            th:=th-t0;

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
                 ti0[1]:=abp;
                 ti0[2]:=ev;
                 ti0[3]:=th;
                 ti0[4]:=prd;
                 ti0[5]:=aa;
                 ti0[6]:=bb;
                 ti0[7]:=ff;
                 ti0[8]:=gg;
            for id:=1 to nt do aex^[id]:=ro^[id];
            sort(aex,nt);
            srdro:=aex^[trunc((nt+1)/2)];
            for id:=1 to nt do
            begin
                mpdq.rtti(ti0,et^[id],rog,tetag);
                dro:=ro^[id]-rog;
                raz:=(teta^[id]-tetag);
                ogr2(raz);
                raz:=ro^[id]*raz*pi/180;
                raz:=sqrt(dro*dro+raz*raz);
                aex^[id]:=raz;
            end;
            sort(aex,nt);
            raz:=aex^[trunc((nt+1)/2)];
     if not ((swds='14575-2125') and (comp10='Ba,Bb')) then
     begin


            if ((raz/srdro)<0.25) or
               (swds='01337-1213') or

               (swds='01361-2954') or
               (swds='06482+5542') or
               (swds='08062+0201') or
               ((swds='08122+1739') and (comp10='Ca,Cb')) or
               (swds='12158-2321') or
               (swds='12335+0901') or
               (swds='21186+1134') or
               ((swds='22038+6438') and (comp10='Aa,Ab'))


               then
            begin
                 nt2:=0;
                 sum:=0;
                 sum2:=0;

                 for id:=1 to nt do
                 begin
                      mpdq.rtti(ti0,et^[id],rog,tetag);
                      dro:=ro^[id]-rog;
                      raz:=(teta^[id]-tetag);
                      ogr2(raz);
                      if ((-220<raz) and (raz<-140)) or( (140<raz) and (raz<220)) then
                      begin
                           teta^[id]:=teta^[id]+180;
                           if teta^[id]>360 then teta^[id]:=teta^[id]-360;
                           x^[id]:=ro^[id]*sin(teta^[id]*pi/180);
                           y^[id]:=ro^[id]*cos(teta^[id]*pi/180);
                      end;
                 end;
            end;
     end else
     begin
                 nt2:=0;
                 sum:=0;
                 sum2:=0;

                 for id:=1 to nt do
                 begin
                      mpdq.rtti(ti0,et^[id],rog,tetag);
                      dro:=ro^[id]-rog;
                      raz:=(teta^[id]-tetag);
                      ogr2(raz);
                      if ((-270<raz) and (raz<-90)) or( (90<raz) and (raz<270)) then
                      begin
                           teta^[id]:=teta^[id]+180;
                           if teta^[id]>360 then teta^[id]:=teta^[id]-360;
                           x^[id]:=ro^[id]*sin(teta^[id]*pi/180);
                           y^[id]:=ro^[id]*cos(teta^[id]*pi/180);
                      end;
                 end;

     end;
     end;


end;
procedure chuzhel2;
var abp,ev,th,prd,aa,bb,ff,gg,comm,somm,comb,somb,ci,si,sum,raz,dte,xa,ya
,xg,yg,sum2:extended;
ti0:das;
nt2:integer;
ai,omb,omm:extended;
begin
{
                  ------P-------  ---a------  --i-----  Node----  -----T-------  -e------  omega---  G  Reference
STF2294              344.837   y    0.732  a  101.60     91.10     1890.086   y  0.935     173.80    4  Sca2015b
                  ------P-------  ---a------  --i-----  Node----  -----T-------  -e------  omega---  G  Reference
J   838              239.84    y    7.451  a   77.6     103.0      1967.51    y  0.784     258.1     5  Ole2002b

HU  628BC             39.4     y    0.70   a  111.6      24.5      1999.1     y  0.35      354.2     3  DRs2012
               +/-     0.1          0.01        1.1       1.3         0.7        0.02        4.6


}
     prd:=39.4;
     abp:=0.70;
     ev:=0.35;
     ai:=111.6;
     omb:=24.5;
     omm:=354.2 ;
     th:=1999.1-t0;
     if  ((swds='06003-3102') and (pos('CE',comp10)=1)) then
     begin
{
                  ------P-------  ---a------  --i-----  Node----  -----T-------  -e------  omega---  G  Reference
TOK   9CE             22.95    y    0.441  a   97.8     146.6      2016.50    y  0.186     193.      3  Tok2014a
               +/-     3.21         0.031       1.0       0.9         0.70       0.094      19.}

             prd:=22.95;
             abp:=0.441;
             ai:=97.8;
             omb:=146.6;
             th:=2016.50-t0;
             ev:=0.186;
             omm:=193 ;
     end;

     if ((swds='08122+1739') and (pos('Ca,Cb',comp10)=1)) then
     begin
{
                  ------P-------  ---a------  --i-----  Node----  -----T-------  -e------  omega---  G  Reference
HUT   1Ca,Cb          19.05    y    0.3728 a  149.2     262.         1983.05  y  0.088      41.      4  RAO2015
               +/-     0.46         0.0068      7.2      12.          0.32       0.013      23.}
             prd:=19.05;
             abp:=0.3728;
             ai:=149.2;
             omb:=262;
             th:=1983.05-t0;
             ev:=0.088;
             omm:=41 ;
     end;



     writeln(rsv);
     writeln(rsv);
     writeln(rsv);
     writeln(rsv);
     writeln(rsv,'╫єцшх ¤ыхьхэЄ√:');
     {prd:=344.837;
     abp:=0.732;
     ev:=0.935;
     ai:=101.60;
     omb:=91.10;
     omm:=173.80;
     th:=1890.086-t0;}

     {prd:=239.84;
     abp:=7.451;
     ev:=0.784;
     ai:=77.6;
     omb:=103.0;
     omm:=258.1 ;
     th:=1967.51-t0;}


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
                 ti0[1]:=abp;
                 ti0[2]:=ev;
                 ti0[3]:=th;
                 ti0[4]:=prd;
                 ti0[5]:=aa;
                 ti0[6]:=bb;
                 ti0[7]:=ff;
                 ti0[8]:=gg;
      nt2:=0;
      sum:=0;
      sum2:=0;

      for id:=1 to nt do
      begin
                mpdq.rtti(ti0,et^[id],rog,tetag);
                dro:=ro^[id]-rog;
                raz:=(teta^[id]-tetag);
                if (raz<-140) or (raz>140) then
                begin
                     teta^[id]:=teta^[id]+180;
                     if teta^[id]>360 then teta^[id]:=teta^[id]-360;
                     x^[id]:=ro^[id]*sin(teta^[id]*pi/180);
                     y^[id]:=ro^[id]*cos(teta^[id]*pi/180);
                end;
                ogr2(raz);
                dte:=raz*pi*ro^[id]/180;
                sum2:=sum2+dro*dro+dte*dte;
           if gbop^[id]=0 then
           begin
                nt2:=nt2+1;
                sum:=sum+dro*dro+dte*dte;
           end;
      end;
      sum:=0.5*sum/nt2;
      sum:=sqrt(sum);
      sum2:=0.5*sum2/nt;
      sum2:=sqrt(sum2);

      writeln(rsv,'╬°.хф.тхё=',sum:22:15);
      writeln(rsv,'╬°.хф.тхё ');
      writeln(rsv,'bez filtr=',sum2:22:15);


      for id:=1 to nt do
      begin

           mpdq.rtti(ti0,et^[id],rog,tetag);
           if gbop^[id]=0 then
           begin

           dro:=ro^[id]-rog;
           dte:=(teta^[id]-tetag);
           ogr2(dte);
           xa:=ro^[id]*sin(teta^[id]*pi/180);
           ya:=ro^[id]*cos(teta^[id]*pi/180);
           xg:=rog*sin(tetag*pi/180);
           yg:=rog*cos(tetag*pi/180);
           ogr2(DTE);

           writeln(rsv,et^[id]+t0:8:3,' ',gbop^[id]:1,' ',gbop3^[id]:1,' '
           ,ro^[id]:8:3,' ',teta^[id]:8:3,' ',rog:8:3,' ',tetag:8:3,' '
           ,dro:8:3,' ',dte:8:3,' ',xa:8:3,' ',ya:8:3);
           end;

      end;
      writeln(rsv);
      for id:=1 to nt do
      begin
           if gbop^[id]<>0 then
           begin
                mpdq.rtti(ti0,et^[id],rog,tetag);

                dro:=ro^[id]-rog;
                dte:=(teta^[id]-tetag);
                ogr2(dte);
                xa:=ro^[id]*sin(teta^[id]*pi/180);
                ya:=ro^[id]*cos(teta^[id]*pi/180);
                xg:=rog*sin(tetag*pi/180);
                yg:=rog*cos(tetag*pi/180);
                ogr2(DTE);

                writeln(rsv,et^[id]+t0:8:3,' ',gbop^[id]:1,' ',ro^[id]:8:3,' ',teta^[id]:8:3,' ',rog:8:3,' ',tetag:8:3,' '
               ,dro:8:3,' ',dte:8:3,' ',xa:8:3,' ',ya:8:3);

           end;

      end;
      flush(rsv);




end;

procedure chuzhel;
var abp,ev,th,prd,aa,bb,ff,gg,comm,somm,comb,somb,ci,si,sum,raz,dte,x,y
,xg,yg,sum2:extended;
ti0:das;
nt2:integer;

begin
{
                  ------P-------  ---a------  --i-----  Node----  -----T-------  -e------  omega---  G  Reference
STF2294              344.837   y    0.732  a  101.60     91.10     1890.086   y  0.935     173.80    4  Sca2015b
                  ------P-------  ---a------  --i-----  Node----  -----T-------  -e------  omega---  G  Reference
J   838              239.84    y    7.451  a   77.6     103.0      1967.51    y  0.784     258.1     5  Ole2002b

}

     writeln(rsv);
     writeln(rsv);
     writeln(rsv);
     writeln(rsv);
     writeln(rsv,'╫єцшх ¤ыхьхэЄ√:');
     {prd:=344.837;
     abp:=0.732;
     ev:=0.935;
     ai:=101.60;
     omb:=91.10;
     omm:=173.80;
     th:=1890.086-t0;}
     prd:=239.84;
     abp:=7.451;
     ev:=0.784;
     ai:=77.6;
     omb:=103.0;
     omm:=258.1 ;
     th:=1967.51-t0;


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
                 ti0[1]:=abp;
                 ti0[2]:=ev;
                 ti0[3]:=th;
                 ti0[4]:=prd;
                 ti0[5]:=aa;
                 ti0[6]:=bb;
                 ti0[7]:=ff;
                 ti0[8]:=gg;
      nt2:=0;
      sum:=0;
      sum2:=0;

      for id:=1 to nt do
      begin
                mpdq.rtti(ti0,et^[id],rog,tetag);
                dro:=ro^[id]-rog;
                raz:=(teta^[id]-tetag);
                ogr2(raz);
                dte:=raz*pi*ro^[id]/180;
                sum2:=sum2+dro*dro+dte*dte;
           if gbop^[id]=0 then
           begin
                nt2:=nt2+1;
                sum:=sum+dro*dro+dte*dte;
           end;
      end;
      sum:=0.5*sum/nt2;
      sum:=sqrt(sum);
      sum2:=0.5*sum2/nt;
      sum2:=sqrt(sum2);

      writeln(rsv,'╬°.хф.тхё=',sum:22:15);
      writeln(rsv,'╬°.хф.тхё ');
      writeln(rsv,'bez filtr=',sum2:22:15);


      for id:=1 to nt do
      begin

           mpdq.rtti(ti0,et^[id],rog,tetag);
           if gbop^[id]=0 then
           begin

           dro:=ro^[id]-rog;
           dte:=(teta^[id]-tetag);
           ogr2(dte);
           x:=ro^[id]*sin(teta^[id]*pi/180);
           y:=ro^[id]*cos(teta^[id]*pi/180);
           xg:=rog*sin(tetag*pi/180);
           yg:=rog*cos(tetag*pi/180);
           ogr2(DTE);

           writeln(rsv,et^[id]+t0:8:3,' ',gbop^[id]:1,' ',gbop3^[id]:1,' '
           ,ro^[id]:8:3,' ',teta^[id]:8:3,' ',rog:8:3,' ',tetag:8:3,' '
           ,dro:8:3,' ',dte:8:3,' ',x:8:3,' ',y:8:3);
           end;

      end;
      writeln(rsv);
      for id:=1 to nt do
      begin
           if gbop^[id]<>0 then
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

                writeln(rsv,et^[id]+t0:8:3,' ',gbop^[id]:1,' ',ro^[id]:8:3,' ',teta^[id]:8:3,' ',rog:8:3,' ',tetag:8:3,' '
               ,dro:8:3,' ',dte:8:3,' ',x:8:3,' ',y:8:3);

           end;

      end;
      flush(rsv);




end;

(*procedure odflM1(pst:string);
const t3=3;
var id,jd,kd,kt,i,mi,mi2,mi3,mi4,nt2,j,t10:integer;
min,max:array[1..t3] of extended;
bmx:boolean;
mix,mix2:array[1..t3] of integer;
rop,st0,mnt,mxt,mxr,raz,raz2,xc,yc,mnt2,svvr,th2:extended;
ti0:das;
t2,diro0,diteta0,drot,dtetat,dr,dvr,dmassa,oev,mint,mprd,mev,mth:extended;
rs,rs2:text;
aa2,bb2,ff2,gg2,comb,somb,comm,somm,ci,f0,tt100:extended;
begin
      prf;
      isuge:=false;
      tfl:=pst;
      niskl:=0;
      ClrScr;
      ishdata(pst);
      if (nt>10)then
      begin
           {clan(m1,m2,teta,x,y,ro,nnt,ain,aex,aex2,gbop,gbop2,gbop3,aex3);}
           for id:=1 to nt do
           begin
                 gbop^[id]:=0;
                 gbop3^[id]:=0;
           end;
           for jd:=1 to t3 do
           begin
                 max[jd]:=-1e39;
                 min[jd]:=1e39;
                 for id:=1 to nt do
                 begin
                      if ro^[id]>max[jd] then
                      begin
                           bmx:=true;
                           for kd:=1 to jd-1 do
                           begin
                                bmx:=bmx and (ro^[id]<max[kd]);
                           end;
                           if bmx then
                           begin
                                max[jd]:=ro^[id];
                                mix[jd]:=id;
                           end;
                      end;
                      if ro^[id]<min[jd] then
                      begin
                           bmx:=true;
                           for kd:=1 to jd-1 do
                           begin
                                bmx:=bmx and (ro^[id]>min[kd]);
                           end;
                           if bmx then
                           begin
                                min[jd]:=ro^[id];
                                mix2[jd]:=id;
                           end;
                      end;
                 end;
           end;
           rop:=(max[t3]+min[t3])*0.5;
           {rop:=min[t3]+(max[t3]-min[t3])/3;}
           for id:=1 to nt do
           begin
                if ro^[id]<rop then
                begin
                      gbop^[id]:=8;
                end;
           end;
           telclan3(et,nnt,ain,aex,aex2,gbop,gbop2,gbop3,5);
           {vydgbop;}
                     st0:=t0;
                     kt:=0;
                     sum:=0;
                     for i:=1 to nnt do
                     begin
                          if gbop^[i]=0 then
                          begin
                               sum:=sum+et^[i]+t0;
                               kt:=kt+1;
                          end;
                     end;
                     sum:=sum/kt;
                     t0:=sum;
                     {t0:=1949.42;}
                     {t0:=1991.25;}
                     for i:=1 to nnt do
                     begin
                          et^[i]:=et^[i]+st0-t0;
                     end;
                nachprib;


           {ciriteracii(true);}
           {mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);}

           mbcir:=true;
           vydres4;
           mbcir:=false;
           mnt:=1e39;
           mxt:=-1e39;
           for id:=1 to nt do
           begin
                if gbop^[id]=0 then
                begin
                  if et^[id]>mxt then mxt:=et^[id];
                  if et^[id]<mnt then mnt:=et^[id];
                end;
           end;
           mi:=1;
           mxr:=-1e39;
           for id:=1 to nt do
           begin
                 if not( (mnt<=et^[id]) and (et^[id]<=mxt)) then
                 begin
                      raz:=abs(mnt-et^[id]);
                      raz2:=abs(mxt-et^[id]);
                      if raz2<raz then raz:=raz2;
                      if raz>mxr then
                      begin
                           mxr:=raz;
                           mi:=id;
                      end;
                 end;
           end;
           mi:=mi-2;
           if mi<1 then mi:=1;
           if mi>nt-5 then mi:=nt-5;
           for id:=0 to 4 do
           begin
                aex^[id+1]:=x^[id+mi];
           end;
           sort(aex,5);
           XC:=aex^[3];
           for id:=0 to 4 do
           begin
                aex^[id+1]:=y^[id+mi];
           end;
           sort(aex,5);
           yc:=aex^[3];
           mnt2:=1e39;
           for id:=0 to 4 do
           begin
                raz:=sqrt(sqr(xc-x^[id+mi])+sqr(yc-y^[id+mi]));
                if raz<mnt2 then
                begin
                     mnt2:=raz;
                     mi2:=id+mi;
                end;
           end;
           kd:=0;
           for id:=1 to nt do
           begin
                if gbop^[id]=0 then
                begin
                     kd:=kd+1;
                     ain^[kd]:=id;
                end;
           end;
           mi3:=ain^[trunc(kd*0.25)];
           mi4:=ain^[trunc(kd*0.75)];
                     {st0:=t0;
                     t0:=et^[mix2[t3]]+st0;
                     for i:=1 to nnt do
                     begin
                          et^[i]:=et^[i]+st0-t0;
                     end;}

           {bro0:=false;
           bteta0:=false;
           btetat:=false;
           brot:=false;}

           t2:=et^[mix2[t3]];
           DrgEphPvd(iro0,iteta0,rot,tetat,r,vr,massa
                    ,t2,diro0,diteta0,drot,dtetat,dr,dvr,dmassa);
           iro0:=diro0;
           iteta0:=diteta0;
           rot:=drot;
           tetat:=dtetat;
           r:=dr;
           vr:=dvr;

                     st0:=t0;
                     t0:=t2+t0;
                     for i:=1 to nnt do
                     begin
                          et^[i]:=et^[i]+st0-t0;
                     end;
           mbcir:=true;
           vydres4;
           mbcir:=false;



           for id:=1 to nt do
           begin
                gbop^[id]:=8;
           end;
           gbop^[mix2[t3]]:=0;
           gbop^[mi2]:=0;
           gbop^[mi3]:=0;
           gbop^[mi4]:=0;
           for id:=1 to nt do
           begin
                gbop^[id]:=0;
           end;


           {gbop^[4]:=0;
           gbop^[14]:=0;}
           uabp:=true;
           uev:=true;
           uth:=false;
           uprd:=true;
           uai:=true;
           uomm:=true;
           uomb:=true;

           bvr:=true;
           bmassa:=true;

           l:=1e-9;
           {riteracii(false);}

           {uth:=true;}
           mbcir:=false;
           r:=r*2;
           if mbcir then
           begin
                mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);
           end else
           begin
                mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);
           end;
                 aa:=ti0[5];
                 bb:=ti0[6];
                 ff:=ti0[7];
                 gg:=ti0[8];

                 {ti[9]:=omm;
                 ti[10]:=omb;
                 ti[11]:=ai;
                 ti[12]:=th;}



            GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
            abp:=ti0[1];
            ev:=ti0[2];
            th:=ti0[3];
            prd:=ti0[4];
            t2:=0;
            GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,t2,diro0,diteta0,drot,dtetat,dr,dvr,dmassa);
            dvr:=GetZnVr(diro0,diteta0,drot,dtetat,dr,dvr,dmassa,pia);

            assign(rs2,'zavev.res');
            rewrite(rs2);
            uth:=true;
            mint:=1e39;
            ev:=j*0.05;
            if ev>=1 then ev:=0.999;
            tt100:=100;
            if (pst12='S01350S2.955') or (pst12='S17066N0.039') then tt100:=10;
            for i:=1 to 500 do
            begin
                 ev:=random;
                 {ev:=0.784;}
                 prd:=random*tt100;
                 th:=random*prd-prd*0.5;
                 GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
                 roo^[i]:=oev;
                 aex^[i]:=ev;
                 aex2^[i]:=prd;
                 aex3^[i]:=th;

                 if (oev<mint)  then
                 begin

                      mint:=oev;
                      mprd:=prd;
                      mev:=ev;
                      mth:=th;
                 end;
           end;
           sort6(roo,aex,aex2,aex3,i);
           assign(rs,'tsmnt.res');
           rewrite(rs);
           mint:=1e39;
           t10:=10;
           if pst12='S08122N1.739' then t10:=1;

           for i:=1 to t10 do
           begin
                ev:=aex^[i];
                prd:=aex2^[i];
                th:=aex3^[i];
                elemiter2(false);
                {readln;}
                writeln(rs,i:3,' ',gf0:8:5,' ',ev:8:5,' ',prd:15:5,' ',th+t0:15:5);
                 if (gf0<mint)  then
                 begin
                      mint:=gf0;
                      mprd:=prd;
                      mev:=ev;
                      mth:=th;
                 end;

           end;
           prd:=mprd;
           ev:=mev;
           th:=mth;


           GetPrmEl(mev,mth,mprd,aa,bb,ff,gg,oev);
           GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
           GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);
           vr:=GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia);
           f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);

           {vr:=-vr;}
           mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);
           {if abs(ti0[2]-ev)>0.001 then vr:=-vr;}
           ClrScr;

           close(rs);
           {halt(0);}
           vydres4;


           bvr:=true;
           bmassa:=true;
           {l:=1e-9;
           riteracii(false);
           vydres4;
           svvr:=vr;}

           {halt(0);}

           for id:=1 to nt do
           begin
                gbop^[id]:=0;
           end;
           {bro0:=true;
           bteta0:=true;
           btetat:=true;
           brot:=true;}
           {koi:=false;}
           {t50000:=1000;}
           {nt2:=nt;
           for id:=1 to nt do
           begin
                if (et^[id]+t0)>2002 then
                begin
                     for j:=1 to 3 do
                     begin
                          nt2:=nt2+1;
                          et^[nt2]:=et^[id];
                          ro^[nt2]:=ro^[id];
                          teta^[nt2]:=teta^[id];
                          x^[nt2]:=x^[id];
                          y^[nt2]:=y^[id];
                          gbop^[nt2]:=0;
                          gbop3^[nt2]:=0;
                     end;
                end;
           end;
           nt:=nt2;
           nnt:=nt;}


           {riteracii(FALSE);}
           elemiter2(true);
           GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
           GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
           GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);
           vr:=GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia);

           {VR:=-vr;}

           th2:=th+t0;
           while th2<(2000-prd) do th2:=th2+prd;
           while th2>(2000) do th2:=th2-prd;
           writeln(rs2,ev:10:5,' ',gf0:10:5,' ',th2:13:5,' ',prd:13:5);
           flush(rs2);
           {elemiter(true);}
           {GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);}
           {vr:=-vr;}

           vydres4;
           {chuzhel;}
           i:=0;
           for id:=1 to nt do
           begin
                if gbop^[id]<>0 then i:=i+1;
           end;

           writeln(rstr,swds,' ',comp10,' ',abp:20:15,' ',ev:20:15,' ',th+t0:20:15,' ',
                         prd:27:15,' ',ai:20:15,' ',omm:20:15,' ',omb:20:15,' ',
                         t0:20:15,' ',
                         iro0:20:15,' ',iteta0:20:15,' ',rot:20:15,' ',
                         tetat:20:15,' ',r:20:15,' ',vr:20:15,' ',massa:20:15,' ',gf0:20:15,' ',nt:4,' ',i:4);
           close(rs2);
           {koi:=false;
           t50000:=300;}
           {riteracii(true);}
           {for id:=1 to nt do
           begin
                if abs(et^[id]+t0-1998.673)<0.01 then
                begin
                     teta^[id]:=teta^[id]-180;
                end;

           end;


           vydres4;
           chuzhel;}


           {oprza:=not(oprza);
           vr:=-vr;
           vydres4;}


           {massa:=massa*2;
           bmassa:=false;
           riteracii(false);
           vydres4;

           bmassa:=true;
           riteracii(false);
           vydres4;}




           kd:=0;


      end;



end;
procedure odflM3(pst:string);
begin

      prf;
      isuge:=false;
      tfl:=pst;
      niskl:=0;
      ClrScr;
      ishdata(pst);
      odflM3d(pst);
end;*)
procedure odflM4(pst:string);
var rs2:text;
begin

      {prf;}
      isuge:=false;
      tfl:=pst;
      niskl:=0;
      ClrScr;
      ishdata2(pst);
      if nt>0 then ish26r;
      if nt>t1330 then
      begin
           ishdata3(pst);
           if nt>0 then ish26r3;

      end;
      {getpvd;}
      if bVneIskl then
      begin
           GetVneIskl;
           if not(bVneIskl2) and (nt>6) then
           begin
                assign(rs2,'errvneis.are');
                mappend(rs2);
                writeln(rs2,'error fl=',pst,' swds=',swds,' comp=',comp10);
                close(rs2);
           end;
      end;

      if nt>3 then
      begin
           getves;
           chuzhel3;
           if opit>0.00001 then
           begin
                odflM5d(pst);
           end else
           begin
                odflM4d(pst);
           end;
      end;
end;

(*procedure odflM2(pst:string);
const t3=3;
var id,jd,kd,kt,i,mi,mi2,mi3,mi4,nt2,j,t10:integer;
min,max:array[1..t3] of extended;
bmx:boolean;
mix,mix2:array[1..t3] of integer;
rop,st0,mnt,mxt,mxr,raz,raz2,xc,yc,mnt2,svvr,th2:extended;
ti0:das;
t2,diro0,diteta0,drot,dtetat,dr,dvr,dmassa,oev,mint,mprd,mev,mth:extended;
rs,rs2:text;
aa2,bb2,ff2,gg2,comb,somb,comm,somm,ci,f0,tt100,f1:extended;
bcm:boolean;
cai,tai,clomb,tomb:extended;
begin
      prf;
      isuge:=false;
      tfl:=pst;
      niskl:=0;
      ClrScr;
      ishdata(pst);
      mi:=-1;
      FOR I:=1 TO 15 do
      begin
           if swds=awdst[i] then mi:=i;
      end;
      pst2:=copy(comp10,1,5);
      bcm:=false;
      if mi>0 then
      begin
            bcm:=(pst2=acompt[mi,1]) or (pst2=acompt[mi,2]);
      end;
      if {(nt>10) and} (mi>0) and bcm then
      begin
           {clan(m1,m2,teta,x,y,ro,nnt,ain,aex,aex2,gbop,gbop2,gbop3,aex3);}
           for id:=1 to nt do
           begin
                 gbop^[id]:=0;
                 gbop3^[id]:=0;
           end;
           uabp:=true;
           uev:=true;
           uth:=false;
           uprd:=true;
           uai:=true;
           uomm:=true;
           uomb:=true;

           bvr:=true;
           bmassa:=true;

            assign(rs2,'zavev.res');
            rewrite(rs2);
            uth:=true;
            mint:=1e39;
            ev:=j*0.05;
            if ev>=1 then ev:=0.999;
            tt100:=2000;

            {if (pst12='S01350S2.955') or (pst12='S17066N0.039') then tt100:=10;}
            for i:=1 to 500 do
            begin
                 ev:=random;
                 {ev:=0.784;}
                 prd:=random*tt100;
                 th:=random*prd-prd*0.5;
                 GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
                 roo^[i]:=oev;
                 aex^[i]:=ev;
                 aex2^[i]:=prd;
                 aex3^[i]:=th;

                 if (oev<mint)  then
                 begin

                      mint:=oev;
                      mprd:=prd;
                      mev:=ev;
                      mth:=th;
                 end;
           end;
           sort6(roo,aex,aex2,aex3,i);
           assign(rs,'tsmnt.res');
           mappend(rs);
           mint:=1e39;
           t10:=10;
           {if pst12='S08122N1.739' then t10:=1;}
           {t10:=1;}

           for i:=1 to t10 do
           begin
                ev:=aex^[i];
                prd:=aex2^[i];
                th:=aex3^[i];
                elemiter2(false);
                {readln;}
                writeln(rs,swds,' ',comp10,' ',i:3,' ',gf0:8:5,' ',ev:8:5,' ',prd:15:5,' ',th+t0:15:5);
                 if (gf0<mint)  then
                 begin
                      mint:=gf0;
                      mprd:=prd;
                      mev:=ev;
                      mth:=th;
                 end;

           end;
           prd:=mprd;
           ev:=mev;
           th:=mth;


           GetPrmEl(mev,mth,mprd,aa,bb,ff,gg,oev);
           GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
           GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);
           vr:=GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia);
           f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);

           mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);
           ClrScr;

           close(rs);


           bvr:=true;
           bmassa:=true;
           for id:=1 to nt do
           begin
                gbop^[id]:=0;
           end;
           {riteracii(FALSE);}
           elemiter2(true);
           GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
           GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
           GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);
           vr:=GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia);

           th2:=th+t0;
           while th2<(2000-prd) do th2:=th2+prd;
           while th2>(2000) do th2:=th2-prd;
           writeln(rs2,ev:10:5,' ',gf0:10:5,' ',th2:13:5,' ',prd:13:5);
           flush(rs2);
           vydnev(1);
           vydres4;
           {chuzhel;}
           i:=0;
           for id:=1 to nt do
           begin
                if gbop^[id]<>0 then i:=i+1;
           end;

           writeln(rstr,swds,' ',comp10,' ',abp:20:15,' ',ev:20:15,' ',th+t0:20:15,' ',
                         prd:27:15,' ',ai:20:15,' ',omm:20:15,' ',omb:20:15,' ',
                         t0:20:15,' ',
                         iro0:20:15,' ',iteta0:20:15,' ',rot:20:15,' ',
                         tetat:20:15,' ',r:20:15,' ',vr:20:15,' ',massa:20:15,' ',gf0:20:15,' ',nt:4,' ',i:4);
           flush(rstr);
           close(rs2);

           mi:=-1;
           FOR I:=1 TO 15 do
           begin
                if swds=atswds[i] then mi:=i;
           end;

           {f0:=fnev3(abp,ev,th,prd,ai,omm+pi,omb+pi);
           f1:=fnev3(abp,ev,th,prd,ai,omm+pi,omb+pi);}
           {elemiter(false);}
           uomb:=false;
           uai:=false;

           h:=1;
           clomb:=atomb[mi]*180/pi;
           tomb:=omb*180/pi;
           raz:=abs(clomb-tomb);
           raz2:=abs(clomb-180-tomb);
           ogr2(raz);
           ogr2(raz2);
           raz:=abs(raz);
           raz2:=abs(raz2);
           if raz2<raz then
           begin
                clomb:=clomb-180;
           end;
           if clomb<0 then clomb:=clomb+360;

           {tomb:=clomb;
           omb:=tomb*pi/180;
           cai:=atai[mi]*180/pi;
           tai:=ai*180/pi;
           if ((cai>90) and (tai<90)) or ((cai<90) and (tai>90)) then
           begin
                cai:=180-cai;
           end;
           ai:=cai*pi/180;
           elemiter(false);}


           if clomb<tomb then h:=-1;
           if abs(clomb-tomb)>180 then h:=-h;
           raz:=clomb-tomb;
           ogr2(raz);
           while abs(raz)>2 do
           begin
                tomb:=tomb+h;
                if tomb<0 then tomb:=tomb+360;
                if tomb>360 then tomb:=tomb-360;
                raz:=clomb-tomb;
                ogr2(raz);
                omb:=tomb*pi/180;
                elemiter(false);
                I:=0;
           end;
           tomb:=clomb;
           omb:=tomb*pi/180;

           elemiter(false);


           uai:=false;


           cai:=atai[mi]*180/pi;
           tai:=ai*180/pi;
           if ((cai>90) and (tai<90)) or ((cai<90) and (tai>90)) then
           begin
                cai:=180-cai;
           end;


           h:=1;
           if cai<tai then h:=-1;
           while abs(cai-tai)>2 do
           begin
                tai:=tai+h;
                ai:=tai*pi/180;
                elemiter(false);
           end;
           tai:=cai;
           ai:=tai*pi/180;

           elemiter(false);
           GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);
           vr:=GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia);
           i:=0;
           for id:=1 to nt do
           begin
                if gbop^[id]<>0 then i:=i+1;
           end;
           vydres4;


           writeln(rstr2,swds,' ',comp10,' ',abp:20:15,' ',ev:20:15,' ',th+t0:20:15,' ',

                         prd:27:15,' ',ai:20:15,' ',omm:20:15,' ',omb:20:15,' ',
                         t0:20:15,' ',
                         iro0:20:15,' ',iteta0:20:15,' ',rot:20:15,' ',
                         tetat:20:15,' ',r:20:15,' ',vr:20:15,' ',massa:20:15,' ',gf0:20:15,' ',nt:4,' ',i:4);
           flush(rstr2);



           kd:=0;


      end;



end;*)


procedure odfl(pst:string);
var {sros:extended;}
ib,is:integer;
siro0,siteta0,srot,stetat,sr,svmassa,svr,dgf0,st0,sum,kt:extended;
ti0:das;


begin
      prf;
      isuge:=false;
      tfl:=pst;
      niskl:=0;
      ClrScr;
      ishdata(pst);
      {regular;}
      sros:=0;
      if nt>1 then
      begin
            sros:=sro/nt;
      end;
      writeln(rsv);
      writeln(rsv);
      writeln(rsv);
      writeln(rsv,pst,' ',swds,' ',comp10,' ',nt:4,' ',sros:10:5);
      flush(rsv);

      {if sros<0.7 then bprd:=false;
      if not(bprd) then
      begin
           exit;
      end;
      if sorb1>5 then
      begin
           if sorb2/sorb1<0.3 then
           begin
                exit;
           end;
      end;}

      if (nt>10) {and (nt<10000)} then
      begin
           clan(m1,m2,teta,x,y,ro,nnt,ain,aex,aex2,gbop,gbop2,gbop3,aex3);
           nnt2:=0;
           sros:=0;
           for ib:=1 to nnt do
           begin
               if gbop^[ib]=0 then
               begin
                    nnt2:=nnt2+1;
                    sros:=sros+ro^[ib];
               end;
           end;
           if nnt2>0 then
           begin
                sros:=sros/nnt2;
           end;
           {writeln('iskl=',nnt-nnt2,' nnt=',nnt);}
           {halt(0);}

           {writeln('clan vse');
           readln;}
           {filtrpol(pst);}
           filtrpol(pst,nnt,gbop,x,y,ro,ain,mxk2);
           { if nnt2>10 then}
           begin
                p1:=pos('\AS',pst);
                pst2:=copy(pst,p1+3,length(pst)-p1-2-4);
                write(rso,pst2,' ',swds,' ',comp10,' ',nt:4,' ',sros:10:5,' ');
                if r26 then
                begin
                     write(rso,' 1');
                end else
                begin
                    write(rso,' 0');
                end;
                write(rso,sorb1:4,' ',sorb2:4,' ',sorb3:4,' ',sorb4:4,' ',sorb5:4,' ');
                flush(rso);

                     st0:=t0;
                     kt:=0;
                     sum:=0;
                     for i:=1 to nnt do
                     begin
                          if gbop^[i]=0 then
                          begin
                               sum:=sum+et^[i]+t0;
                               kt:=kt+1;
                          end;
                     end;
                     sum:=sum/kt;
                     t0:=sum;
                     {t0:=1991.25;}
                     for i:=1 to nnt do
                     begin
                          et^[i]:=et^[i]+st0-t0;
                     end;
                nachprib;
                {halt(0);}
                writeln(rsv,'═рўры№эюх яЁшсышцхэшх:');
                if posf then
                begin
                     write(rso,' 9');
                end else
                begin
                    write(rso,' 0');
                end;

                vydres4;
                ClrScr;
                {halt(0);}

                {if itrf<3000 then}
                {if false then}
                writeln('filtr pol mxk2=',mxk2);
                readln;


                if mxk2=0 then
                begin
                     ciriteracii2(true);
                     {gbop^[7]:=1;}
                     if not(bzc) then  ciriteracii(true);

                     {ciriteracii(false);}
                     {readln;}
                     vydnev(1);
                     mbcir:=true;
                     vydres4;
                     mbcir:=false;


                     if (r<999.999) and (uci) then
                     begin
                          writeln(rsv);
                          writeln(rsv);
                          writeln(rsv);
                          bmassa:=true;
                          bvr:=true;

                          ItrIskl1;
                          ItrIskl2;
                          for ib:=1 to nnt do
                          begin
                               if (et^[ib]+t0)=1841.800 then gbop^[ib]:=3;
                          end;

                          riteracii(true);
                          vydnev(2);

                          vydres4;
                     end;
                end

                else
                begin
                     {writeln(rsv,'Ne soshlis iteracii v krugovom dvigrnii');}
                     t50000:=50000;
                     predv(false);

                     for is:=1 to nnt do gbop2^[is]:=gbop^[is];
                     {bvr:=true;}

                     {riteracii(true);}
                     {gbop^[1]:=0;}

                     siro0:=iro0;
                     siteta0:=iteta0;
                     srot:=rot;
                     stetat:=tetat;
                     sr:=r;
                     svmassa:=massa;
                     svr:=vr;



                     ciriteracii(true);
                     dgf0:=gf0;
                     {bmassa:=true;}
                     {writeln('poln itr');
                     writeln('poln itr');
                     writeln('poln itr');
                     readln;}
                     {iro0:=siro0;
                     iteta0:=siteta0;
                     rot:=srot;
                     tetat:=stetat;
                     vr:=vr*0.9;
                     massa:=massa*0.9;}

                     {sr:=r;
                     svmassa:=massa;}
                     {for i:=1 to nnt do
                     begin
                          if gbop^[i]>7 then gbop^[i]:=0;
                     end;}
                     mbcir:=true;
                     vydres4;
                     mbcir:=false;



                     bmassa:=true;
                     bvr:=true;

                     ItrIsklB7;
                     ItrIskl1;
                     ItrIskl2;

                     {elemiter(true);}
                     {ClrScr;}
                     riteracii(true);
                     vydres4;


                     (*iro0:=siro0;
                     iteta0:=siteta0;
                     rot:=srot;
                     tetat:=stetat;
                     r:=sr;
                     massa:=svmassa;
                     vr:=svr;

                     uabp:=true;
                     uev:=true;
                     uth:=true;
                     uprd:=true;
                     uai:=true;
                     uomm:=true;
                     uomb:=true;



                     elemiter(false);*)
                     vydnev(3);

                     {halt(0);}

                     mxk2:=0;


                end;
                writeln(rso);
                {ClrScr;
                iteracii(true);
                writeln(rsv);
                writeln(rsv);
                writeln(rsv);
                writeln(rsv,'╚Єюуютр  юЁсшЄр:');
                vydres4;}
           end;
      end;
end;
function bfln(otm:string):boolean;
var vzv:boolean;
i:integer;
begin
          assign(dt,otm);
          reset(dt);
          vzv:=false;
          while not(eof(dt)) and not(vzv) do
          begin
               pst:='';
               ch:=' ';
               i:=0;
               while (ch<>#10) and (i<256) do
               begin
                    read(dt,ch);
                    pst:=pst+ch;
                    i:=i+1;
               end;
               if
                       (pst[12]='.')
                  and  (pst[23]='.')
                  and  (pst[39]='.')
               then
               begin
                    vzv:=true;
               end;
          end;
          close(dt);
          bfln:=vzv;
end;
function GetIshFile:string;
var de:integer;
    hu:searchrec;
    vzv:string;
    dt:text;


begin
      GetIshFile:='input.txt';
      exit;
      {assign(dt,'odfl.txt');
      reset(dt);
      readln(dt,vzv);
      GetIshFile:=vzv;
      exit;
      close(dt);}


      findfirst('w*.txt',$3F,hu);
      de:=doserror;
      vzv:='';
      if de=0 then
      begin
           {hu.name;}
           if bfln(hu.name) then vzv:=hu.name;

      end;
      while (de=0) and (vzv='') do
      begin
           findnext(hu);
           de:=doserror;
           if bfln(hu.name) then vzv:=hu.name;

     end;
     GetIshFile:=vzv;


end;
function cyf(rf:char):boolean;
var tc:boolean;
begin
     tc:= (rf='1') or (rf='2') or (rf='3') or (rf='4') or (rf='5') or
     (rf='6') or (rf='7') or (rf='8') or (rf='9') or  (rf='0');
     cyf:=tc;
end;
procedure corinpr(fl:string);
var dt1,dt2:text;
i,j:integer;
ep:extended;
begin
     assign(dt1,fl);
     reset(dt1);
     for i:=1 to 12 do readln(dt1);
     for i:=1 to 12 do
     begin
           readln(dt1,ep);
           j:=0;
     end;
     close(dt1);

end;

procedure corinp(fl:string);
var dt,rs:text;
ch,pch:char;
nb,nb2,i:longint;
p:integer;
nfl,rfl,nnfl:string;
f:file;
begin
     p:=pos('.',fl);
     nfl:=fl;
     rfl:='';
     if p>0 then
     begin
          nfl:=copy(fl,1,p-1);
          rfl:=copy(fl,p+1,length(fl)-p);
          i:=0;
     end;
     assign(dt,fl);
     reset(dt);
     nb:=0;
     nb2:=0;

     while not(eof(dt)) do
     begin
           read(dt,ch);
           nb:=nb+1;
           if cyf(ch) then nb2:=nb;

     end;
     close(dt);
     nnfl:=nfl+'.qvk';
     assign(dt,fl);
     reset(dt);

     assign(rs,nnfl);
     rewrite(rs);
     for i:=1 to nb2 do
     begin
          read(dt,ch);
          write(rs,ch);
     end;
     close(dt);
     close(rs);
     assign(dt,nnfl);
     reset(dt);

     nnfl:=nfl+'.qvr';

     assign(rs,nnfl);
     rewrite(rs);
     pch:=' ';
     nb:=0;
     while not(eof(dt)) do
     begin
          read(dt,ch);
          nb:=nb+1;
          if (nb<>1) and (ch=#10)  then
          begin
                if pch<>#13 then
                begin
                     write(rs,#13);
                end;
          end;
          write(rs,ch);
          pch:=ch;

     end;
     close(dt);
     close(rs);
     assign(f,fl);
     erase(f);
     assign(f,nnfl);
     rename(f,fl);
     assign(f,nfl+'.qvk');
     erase(f);

end;


begin
      corinp('input.txt');
      {corinp('w.txt');}
      {halt(0);}
      bklpvd:=false;
      {bVneIskl:=true;}
      bVneIskl:=false;
      model:=false;
      fixoi:=false;
      fixoi2:=false;
      koi:=true;
      bmjd:=tekmjdr;

      {assign(rsfp,'filtrpl.res');
      rewrite(rsfp);}
      {assign(ris,'iskl.res');
      rewrite(ris);
      assign(rsq,'oev.res');
      rewrite(rsq);
      assign(rso,'otbr.res');
      rewrite(rso);}

      oflnm:=GetIshFile;
      mflnm:=oflnm;

     mdf:=0;
     ku:=0;
     ku2:=0;
     iud:=0;
     inud:=0;
     new(ro);
     new(teta);
     new(et);
     new(x);
     new(y);
     new(gbop);
     new(gbop2);
     new(gbop3);
     new(agmi);
     new(aex);
     new(aex2);
     new(aex3);
     new(ain);

     {new(xro0);
     new(xteta0);}

     new(a);
     new(at);
     new(al);
     for i:=1 to mnk.nk2 do new(al^[i]);

      mbcir:=false;
      new(aref);
      {new(aref2);}

      new(sf);
      {new(sf2);}
      new(aves);
      {new(roo);
      new(tetao);}

      {procedure odfl(pst:string);}
      p1:=pos('.',mflnm);
      pst:=copy(mflnm,1,p1)+'nt2';

      {assign(rsv,pst);
      rewrite(rsv);
      close(rsv);}
      p1:=pos('.',mflnm);
      pst:=copy(mflnm,1,p1)+'orb';
      pst:='output.txt';
      assign(rsv,pst);
      rewrite(rsv);
      ncomp:=0;
      repeat
      begin
                pst:=mflnm;
                begin
                     omb:=0;
                     ai:=0;
                     svpst:=pst;
                     ncomp:=ncomp+1;
                     odflM4(pst);

                     {if not(btesorb) then}
                     if nt>3 then
                     begin
                          bmjd:=tekmjdr;

                          fixoi2:=true;
                          pst:=svpst;
                          if not(bitrves) then odflM4d(pst);
                          fixoi2:=false;
                          vydres4{(ev,th,prd,abp,ai,omb,omm)};

                     end;



                end;
            writeln(rsv);
            writeln(rsv);
      end;

      until true {not(ecm)};
      bmjd:=tekmjdr-bmjd;
      bmjd:=bmjd*60*24;
      close(rsv);
      {assign(rs0,oflnm+'j.res');
      rewrite(rs0);
      writeln(rs0,bmjd:10:5);
      close(rs0);
      close(rsq);
      close(ris);
      close(rso);
      close(rsfp);
      close(rstr);
      close(rstr2);
            if not(boshi) and bvne then
            begin

                 close(arsv[1]);
                 close(arsv[2]);
                 close(arsv[3]);
                 close(rsv);
            end;}


      {writeln('ku-ku');
      writeln('ku-ku2');
      writeln('ku-ku3');
      writeln('ku-ku4');
      while i<100000 do
      begin
           i:=i+1;
           j:=i;
      end;}
      {close(rsfp);}





end.