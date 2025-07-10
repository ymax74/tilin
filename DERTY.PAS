unit derty;

interface
         uses sbop,mpdq,ubrk,mnk,mpvd,crt,dos;
         procedure vydres4;
         procedure riteracii(biskl:boolean);
         procedure ciriteracii(biskl:boolean);
         procedure elemiter(biskl:boolean);
         procedure elemiter2(biskl:boolean);
         function fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia:extended):extended;
         function fnev3(abp,ev,th,prd,ai,omm,omb:extended):extended;
         function flam(lam:extended;vgrad:das;iro0,iteta0,rot,tetat,r,vr,massa,pia:extended):extended;
         {procedure cvb(imxr:integer);}
         function GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia:extended):extended;
         procedure  getelem(abp,ev,th,prd,ai,omm,omb,t,pia:extended;var ti:das);
         function fnev4(ev,th,prd:extended):extended;

implementation
(*procedure cvb(imxr:integer);
var dte:extended;
begin
           niskl:=niskl+1;
           aiskl[niskl]:=imxr;
           dte:=(tetao^[imxr]-teta^[imxr]);
           ogr2(dte);
           if not(bvnevyd) then
           begin
                writeln('nabl  N=',imxr:3,' god= ',(et^[imxr]+t0):10:5,' nevazka ro=',roo^[imxr]-ro^[imxr]:10:5,
                ' nevazka te=',dte:10:5);
           end;
           if not(isuge) then writeln(ris,tfl,' ',swds,' ',comp10,' ',nt:4);
           isuge:=true;
           writeln(ris,imxr:3,' ',(et^[imxr]+t0):10:5,' ',roo^[imxr]-ro^[imxr]:10:5,
           ' ',dte:10:5,' ',oevr:10:5,' ',oevt:10:5);
           {ephi^[niskl]:=et^[imxr];
           tetai^[niskl]:=teta^[imxr];}

           {roi^[niskl]:=ro^[imxr];}
           for i:=imxr to nt-1 do
           begin
                et^[i]:=et^[i+1];
                ro^[i]:=ro^[i+1];
                teta^[i]:=teta^[i+1];
                x^[i]:=ro^[i]*sin(teta^[i]*pi/180);
                y^[i]:=ro^[i]*cos(teta^[i]*pi/180);
                aref^[i]:=aref^[i+1];
           end;
           nt:=nt-1;
end;*)

procedure cvb2(imxr:integer;t2:integer);
var dte:extended;
begin
           niskl:=niskl+1;
           aiskl[niskl]:=imxr;
           if not(bvnevyd) then
           begin

                 writeln('N=',imxr:3,' Ep= ',(et^[imxr]+t0):10:5,'  d rho=',gro:10:5,
                 ' d theta=',gte:10:5);
           end;
           {if not(isuge) then writeln(ris,tfl,' ',swds,' ',comp10,' ',nt:4);}
           isuge:=true;
           {writeln(ris,imxr:3,' ',(et^[imxr]+t0):10:5,' ',roo^[imxr]-ro^[imxr]:10:5,
           ' ',dte:10:5,' ',oevr:10:5,' ',oevt:10:5);}
           {ephi^[niskl]:=et^[imxr];
           tetai^[niskl]:=teta^[imxr];
           roi^[niskl]:=ro^[imxr];}
           gbop^[imxr]:=t2;

end;


procedure rvh;
var dt:text;
    idc,jdc,kdc:integer;
begin
     assign(dt,'vh.dt');
     reset(dt);
     readln(dt,iro0);
     readln(dt,iteta0);
     readln(dt,rot);
     readln(dt,tetat);
     readln(dt,r);
     readln(dt,vr);
     readln(dt,massa);
     readln(dt,pia);
     read(dt,jdc);
     niskl:=0;
     for idc:=1 to jdc do
     begin
          read(dt,kdc);
          {cvb(kdc);}
     end;
     readln(dt);
     close(dt);
end;

procedure rvh2;
var dt:text;
    idc,jdc,kdc:integer;
begin
     assign(dt,'vh.dt');
     reset(dt);
     readln(dt,iro0);
     readln(dt,iteta0);
     readln(dt,rot);
     readln(dt,tetat);
     readln(dt,r);
     readln(dt,vr);
     readln(dt,massa);
     readln(dt,pia);
     read(dt,jdc);
     close(dt);
end;

procedure maxiskl(iro0,iteta0,rot,tetat,r,vr,massa,pia:extended;var oevr,oevt:extended;var maxnev:extended;var mi:integer);
var
ti0:das;
dro,dte,sum2,sum:extended;
p1,p2:integer;
nt2:integer;
begin

      mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);
      sum:=0;
      sum2:=0;
      nt2:=0;
      for id:=1 to nt do
      begin
           if gbop^[id]=0 then
           begin
                nt2:=nt2+1;
                mpdq.rtti(ti0,et^[id],rog,tetag);
                dro:=rog-ro^[id];
                aex^[nt2]:=dro;
           end;
      end;
      sort(aex,nt2);
      p1:=trunc(nt*intver);
      p2:=trunc(nt*(1-intver));

      oevr:=(aex^[p2]-aex^[p1])*0.5;
      nt2:=0;
      for id:=1 to nt do
      begin
           if gbop^[id]=0 then
           begin
                nt2:=nt2+1;
                mpdq.rtti(ti0,et^[id],rog,tetag);
                dte:=(tetag-teta^[id]);
                ogr2(dte);
                aex^[nt2]:=dte;
           end;

      end;
      p1:=trunc(nt*intver);
      p2:=trunc(nt*(1-intver));
      sort(aex,nt);
      oevt:=(aex^[p2]-aex^[p1])*0.5;

      maxnev:=0;
      mi:=-1;
      for id:=1 to nt do
      begin
           mpdq.rtti(ti0,et^[id],rog,tetag);
           if gbop^[id]=0 then
           begin
                dro:=abs(rog-ro^[id])/oevr;
                dte:=abs(tetag-teta^[id]){/oevt};
                ogr2(dte);
                dte:=abs(dte)/oevt;
                if id=14 then
                begin
                     id:=14;
                end;
                if dte>dro then dro:=dte;
                if dro>maxnev then
                begin
                     maxnev:=dro;
                     mi:=id;
                end;
           end;
           {roo^[id]:=rog;
           tetao^[id]:=tetag;}

      end;
      gmi:=mi;

      p1:=0;
end;

procedure maxiskl4(ev,th,prd:extended;var oevr,oevt:extended;var maxnev:extended;var mi:integer);
var
ti0:das;
dro,dte,sum2,sum:extended;
p1,p2:integer;
nt2:integer;
aa,bb,ff,gg,oev:extended;
begin
      GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);

      sum:=0;
      sum2:=0;
      nt2:=0;
      for id:=1 to nt do
      begin
           if gbop^[id]=0 then
           begin
                nt2:=nt2+1;
                {mpdq.rtti(ti0,et^[id],rog,tetag);}
                mpdq.rtti3(et^[id],ev,th,prd,aa,bb,ff,gg,rog,tetag);

                dro:=rog-ro^[id];
                dro:=dro*aves^[id];
                aex^[nt2]:=dro;
           end;
      end;
      sort(aex,nt2);

      p1:=trunc(nt*intver);
      p2:=trunc(nt*(1-intver));
      if p1<1 then p1:=1;
      if p1>nt then p1:=nt;
      if p2<1 then p2:=1;
      if p2>nt then p2:=nt;


      oevr:=(aex^[p2]-aex^[p1])*0.5;
      nt2:=0;
      for id:=1 to nt do
      begin
           if gbop^[id]=0 then
           begin
                nt2:=nt2+1;
                {mpdq.rtti(ti0,et^[id],rog,tetag);}
                mpdq.rtti3(et^[id],ev,th,prd,aa,bb,ff,gg,rog,tetag);

                dte:=(tetag-teta^[id])*aves^[id];
                ogr2(dte);
                aex^[nt2]:=dte;
           end;

      end;
      p1:=trunc(nt*intver);
      p2:=trunc(nt*(1-intver));
      sort(aex,nt);
      if p1<1 then p1:=1;
      if p1>nt then p1:=nt;
      if p2<1 then p2:=1;
      if p2>nt then p2:=nt;

      oevt:=(aex^[p2]-aex^[p1])*0.5;

      maxnev:=0;
      mi:=-1;
      for id:=1 to nt do
      begin
           {mpdq.rtti(ti0,et^[id],rog,tetag);}
           mpdq.rtti3(et^[id],ev,th,prd,aa,bb,ff,gg,rog,tetag);

           if gbop^[id]=0 then
           begin
                dro:={aves^[id]*}abs(rog-ro^[id])/oevr;
                dte:={aves^[id]*}abs(tetag-teta^[id]){/oevt};
                ogr2(dte);
                dte:=abs(dte)/oevt;
                if dte>dro then dro:=dte;
                if dro>maxnev then
                begin
                     maxnev:=dro;
                     mi:=id;
                     gte:=abs(tetag-teta^[id]);
                     ogr2(gte);

                     gro:=abs(rog-ro^[id])
                end;
           end;
           {roo^[id]:=rog;
           tetao^[id]:=tetag;}

      end;
      gmi:=mi;

      p1:=0;
end;

procedure wvh;
var dt:text;
    idc:integer;
begin
     {assign(dt,'vh.dt');
     rewrite(dt);
     writeln(dt,iro0);
     writeln(dt,iteta0);
     writeln(dt,rot);
     writeln(dt,tetat);
     writeln(dt,r);
     writeln(dt,vr);
     writeln(dt,massa);
     writeln(dt,pia);
     write(dt,niskl);
     for idc:=1 to niskl do write(dt,' ',aiskl[idc]);
     writeln(dt);
     flush(dt);
     close(dt);}
end;

function flam(lam:extended;vgrad:das;iro0,iteta0,rot,tetat,r,vr,massa,pia:extended):extended;
begin
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


        flam:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);

end;

function fff(lam:extended;vgrad:das;iro0,iteta0,rot,tetat,r,vr,massa,pia:extended):extended;
var f0,f1,eps:extended;
begin
     eps:=0.0001*gllam;

     f0:=flam(lam-eps*0.5,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
     f1:=flam(lam+eps*0.5,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
     fff:=(f1-f0)/eps;
end;
function GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia:extended):extended;
var f0,f1,vvr:extended;
begin
     vvr:=vr;
     f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
     f1:=fnev(iro0,iteta0,rot,tetat,r,-vr,massa,pia);
     if f1<f0 then vvr:=-vr;
     GetZnVr:=vvr;

end;

function fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia:extended):extended;
var
ti0:das;
dro,dte,raz,max:extended;
nt2:integer;
begin
      nevgmi:=-1;
      max:=-1e38;
      if (r<0) or (r>3000) then
      begin
           fnev:=1e38;
           exit;
      end;

      if (massa<0) or (massa>3000) then
      begin
           fnev:=1e38;
           exit;
      end;
      {nt:=1;}
      if mbcir then
      begin
           mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,vr,massa,0,pia,ti0);
           {mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,0,pia,ti0);}
      end else
      begin
           mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);


           if (abs(ti0[2]-1.0)<0.000000001) or (ti0[2]>10) or (ti0[4]<1) then
           begin
                fnev:=1e38;
                exit;

           end;

      end;
      {if ti0[2]>=1 then
      begin
           fnev:=1e39;
      end else}
      begin

      sum:=0;
      nt2:=0;
      {if mdf=1 then
      begin
           dk1:=0
      end;
      if mdf=1 then
      begin
           dk2:=0
      end;}

      for id:=1 to nt do
      begin
           if gbop^[id]=0 then
           begin
                nt2:=nt2+1;
                mpdq.rtti(ti0,et^[id],rog,tetag);
                dro:=rog-ro^[id];
                raz:=(tetag-teta^[id]);
                ogr2(raz);
                dte:=raz*pi*ro^[id]/180;
                {dte:=raz*pi*rog/180;}
                dro:=dro*aves^[id];
                dte:=dte*aves^[id];

                sum:=sum+dro*dro+dte*dte;
                dro:=abs(dro);
                dte:=abs(dte);
                if dro>max then
                begin
                     max:=dro;
                     nevgmi:=id;
                end;
                if dte>max then
                begin
                     max:=dte;
                     nevgmi:=id;

                end;
                {if mdf=1 then
                begin
                     dk1:=dk1+1;
                     aex^[dk1]:=dte;
                end;
                if mdf=2 then
                begin
                     dk2:=dk2+1;
                     aex2^[dk2]:=dte;
                end;}
           end;
      end;
      sum:=0.5*sum/nt2;
      sum:=sqrt(sum);
      fnev:= sum;
      end;
end;

function sign(ve:extended):extended;
begin
     if ve>0 then
     begin
          sign:=1;
     end else
     begin
          sign:=-1;
     end;
end;


procedure vydres4;
var
ti0:das;
dro,dte,x,y,xg,yg,abp,tvr,tmassa,a,tvr2,sum,raz,rob,teb,roe,tee,roc,tec,etc,etb,ete
,sum2,tbk:extended;
rst,rst2:text;
kfl:string;
p1,nt2:integer;
ugeb:boolean;
sumves,kvr,ves:extended;
begin
      {if (nevyh) then
      begin
           writeln(rsv);
           writeln(rsv,'Iteracii ne soshlis');
           exit;

      end;}
      {writeln(rsv,'ku=',ku,' ku2=',ku2);}
      if mbcir then
      begin
           mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,tvr2,tmassa,et^[1],1,ti0);
      end else
      begin
           mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);
      end;
      if not(bklpvd) then
      begin

          if glomb>pi then
          begin
               glomb:=glomb-pi;
               glomm:=glomm-pi;
               if glomm<0 then glomm:=2*pi+glomm;
          end;

            comm:=cos(glomm);
            somm:=sin(glomm);
            comb:=sin(glomb);
            somb:=cos(glomb);
            ci:=cos(glai);
            si:=sin(glai);
            aa:= comm*comb-somm*somb*ci;
            bb:= comm*somb+somm*comb*ci;
            ff:=-somm*comb-comm*somb*ci;
            gg:=-somm*somb+comm*comb*ci;

      ti0[2]:=glev;
      ti0[3]:=glth;
      ti0[4]:=glprd;
      ti0[1]:=glabp;
      ti0[11]:=glai;
      ti0[10]:=glomb;
      ti0[9]:=glomm;
      ti0[5]:=aa;
      ti0[6]:=bb;
      ti0[7]:=ff;
      ti0[8]:=gg;
      end;
      abp:=ti0[1];
      tvr:=vr;
      if mbcir then tvr:=tvr2/pia;
      tmassa:=massa;
      if mbcir then
      begin
           a:=abp/pia;
           tmassa:=a*a*a/(ti0[4]*ti0[4]);
      end;
      writeln(rsv,swds,' ',comp10);

      writeln(rsv,'t0            =',t0:22:15);
      writeln(rsv,'rho0          =',iro0:22:15);
      writeln(rsv,'theta0        =',iteta0:22:15);
      writeln(rsv,'rho t         =',rot:22:15);
      writeln(rsv,'theta t       =',tetat:22:15);
      writeln(rsv,'rad.curvature =',r:22:15);
      write  (rsv,'R V["/year]   =',tvr:22:15);
      if bklpvd then write  (rsv,' ',ovr2:8:4);
      writeln(rsv);
      write(rsv,'mass(pi=1)    =',tmassa:22:15);
      if bklpvd then write  (rsv,' ',omas2:8:4,' pi=',opit:8:4);
      writeln(rsv);


      {writeln(rsv,'pi       =',pia:22:15);}
      writeln(rsv);
      tbk:=180-ti0[11]*180/pi;
      if tbk<0 then tbk:=tbk+360;
      writeln(rsv,'a    =',abp:22:15);
      writeln(rsv,'e    =',ti0[2]:22:15);
      writeln(rsv,'T    =',ti0[3]+t0:22:15);
      writeln(rsv,'P    =',ti0[4]:22:15);

      writeln(rsv,'omega=',ti0[9]*180/pi:22:15);
      writeln(rsv,'Node =',ti0[10]*180/pi:22:15);
      writeln(rsv,'i    =',tbk:22:15);
      writeln(rsv);

      if ti0[2]>10 then exit;
      sum:=0;
      p1:=0;
      if mbcir then
      begin
           kfl:=tfl;
           p1:=pos('.TXT',kfl);
           if p1>8 then
           begin
                kfl:=copy(kfl,1,p1-1);
                kfl:=kfl+'b.res';
                assign(rst,kfl);
                rewrite(rst);
           end;

           kfl:=tfl;
           p1:=pos('.TXT',kfl);
           if p1>8 then
           begin
                kfl:=copy(kfl,1,p1-1);
                kfl:=kfl+'d.res';
                assign(rst2,kfl);
                rewrite(rst2);
           end;
      end;
      sumves:=0;
      nt2:=0;
      for id:=1 to nt do
      begin
           if gbop^[id]=0 then
           begin
                nt2:=nt2+1;
                sumves:=sumves+aves^[id];
           end;
      end;
      kvr:=nt2/sumves;

      nt2:=0;
      sum2:=0;

      for id:=1 to nt do
      begin
                mpdq.rtti(ti0,et^[id],rog,tetag);
                dro:=ro^[id]-rog;
                raz:=(teta^[id]-tetag);
                ogr2(raz);
                dte:=raz*pi*ro^[id]/180;
                sum2:=sum2+dro*dro*aves^[id]*aves^[id]+dte*dte*aves^[id]*aves^[id];
           if gbop^[id]=0 then
           begin
                ves:=kvr*aves^[id];
                nt2:=nt2+1;
                sum:=sum+dro*dro*ves*ves+dte*dte*ves*ves;
           end;
      end;
      sum:=0.5*sum/nt2;
      sum:=sqrt(sum);
      sum2:=0.5*sum2/nt;
      sum2:=sqrt(sum2);

      writeln(rsv,'RMS=',sum:22:15);
      writeln(rsv);
      writeln(rsv,'RMS');
      writeln(rsv,'all obs.=',sum2:22:15);
      writeln(rsv);
      glf0q:=sum;

      {for id:=1 to nt do
      begin
           if gbop^[id]=0 then
           begin
                nt2:=nt2+1;
                mpdq.rtti(ti0,et^[id],rog,tetag);
                dro:=ro^[id]-rog;
                raz:=(teta^[id]-tetag);
                ogr2(raz);
                dte:=raz*pi*ro^[id]/180;
                sum:=sum+dro*dro+dte*dte;
           end;
      end;
      sum:=0.5*sum/nt2;
      sum:=sqrt(sum);
      writeln(rsv,'ќш.ед.вес=',sum:22:15);}
      ugeb:=false;
       writeln(rsv,
       '    Ep       Rho(obs)      Theta(obs)     Rho(c)     Theta(c)       dRho        dTheta          X            Y ');
       {1989.780    0.451   42.700    0.458   36.535   -0.007    6.165    0.306    0.331}


      for id:=1 to nt do
      begin

           mpdq.rtti(ti0,et^[id],rog,tetag);
           if gbop^[id]=0 then
           begin

           if not(ugeb) then
           begin
                rob:=rog;
                teb:=tetag;
                etb:=et^[id];
                ugeb:=true;
           end;
           begin
                roe:=rog;
                tee:=tetag;
                ete:=et^[id];
           end;
           dro:=ro^[id]-rog;
           dte:=(teta^[id]-tetag);
           ogr2(dte);
           x:=ro^[id]*sin(teta^[id]*pi/180);
           y:=ro^[id]*cos(teta^[id]*pi/180);
           xg:=rog*sin(tetag*pi/180);
           yg:=rog*cos(tetag*pi/180);
           ogr2(DTE);

           writeln(rsv,et^[id]+t0:10:5,' '
           ,ro^[id]:12:7,' ',teta^[id]:12:7,' ',rog:12:7,' ',tetag:12:7,' '
           ,dro:12:7,' ',dte:12:7,' ',x:12:7,' ',y:12:7{,' ',aref^[id],' ',aves^[id]:12:7});

           {if mbcir and (p1>8) then
           begin
                writeln(rst,et^[id]+t0:8:3,' ',x:8:3,' ',y:8:3,' ',ro^[id]:8:3,' ',teta^[id]:8:3,' '
                ,dro:8:3,' ',dte:8:3);
                writeln(rst,et^[id]+t0:8:3,' ',xg:8:3,' ',yg:8:3,' ',rog:8:3,' ',tetag:8:3,' '
                ,dro:8:3,' ',dte:8:3);

                writeln(rst2,et^[id]+t0:8:3,' ',gbop^[id]:1,' ',x:8:3,' ',y:8:3,' ',ro^[id]:8:3,' ',teta^[id]:8:3,' '
                ,dro:8:3,' ',dte:8:3);


           end;

           end else
           begin
               if mbcir and (p1>8) then
               begin
                    x:=ro^[id]*sin(teta^[id]*pi/180);
                    y:=ro^[id]*cos(teta^[id]*pi/180);
                    writeln(rst2,et^[id]+t0:8:3,' ',gbop^[id]:1,' ',x:8:3,' ',y:8:3,' ',ro^[id]:8:3,' ',teta^[id]:8:3);
               end;}


           end;

      end;
      writeln(rsv);
      writeln(rsv,'excluded:');
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

                {writeln(rsv,et^[id]+t0:8:3,' ',ro^[id]:8:3,' ',teta^[id]:8:3,' ',rog:8:3,' ',tetag:8:3,' '
               ,dro:8:3,' ',dte:8:3,' ',x:8:3,' ',y:8:3,' ',aref^[id]);}
                writeln(rsv,et^[id]+t0:10:5,' '
               ,ro^[id]:12:7,' ',teta^[id]:12:7,' ',rog:12:7,' ',tetag:12:7,' '
               ,dro:12:7,' ',dte:12:7,' ',x:12:7,' ',y:12:7,' ',aref^[id],' ',aves^[id]:12:7);

           end;

      end;
      etc:=(etb+ete)/2;
      mpdq.rtti(ti0,etc,roc,tec);

           if mbcir and (p1>8) then
           begin
                close(rst);
                close(rst2);
           end;
      if mbcir then
      begin
           if uci then
           begin
                write(rso,' 1');
                write(rso,nt:4,' ',niskl:4);
                write(rso,' ',etb+t0:4:0);
                write(rso,' ',ete+t0:4:0);
                write(rso,' ',rob:7:3);
                write(rso,' ',teb:7:3);
                write(rso,' ',roc:7:3);
                write(rso,' ',tec:7:3);
                write(rso,' ',roe:7:3);
                write(rso,' ',tee:7:3);
                write(rso,' '{ro0      =},iro0:22:15);
                write(rso,' '{teta0    =},iteta0:22:15);
                write(rso,' '{rot      =},rot:22:15);
                write(rso,' '{tetat    =},tetat:22:15);
                write(rso,' '{krivizna =},r:22:15);
                for i:=1 to 6 do write(rso,' ',sgm[i]:22:15);

                write(rso,' '{t0       =},t0:22:15);
                write(rso,' '{vr       =},tvr:22:15);
                write(rso,' '{massa    =},tmassa:22:15);
                write(rso,' '{pi       =},pia:22:15);
                write(rso,' '{b. poluos=},abp:22:15);
                write(rso,' '{eks      =},ti0[2]:22:15);
                write(rso,' '{periastr =},ti0[3]+t0:22:15);
                write(rso,' '{period   =},ti0[4]:22:15);

                write(rso,' '{omm      =},ti0[9]*180/pi:22:15);
                write(rso,' '{omb      =},ti0[10]*180/pi:22:15);
                write(rso,' '{i        =},ti0[11]*180/pi:22:15);

           end else
           begin
                write(rso,' 0');
                write(rso,nt:4,' ',niskl:4);

           end;
           flush(rso);
      end;
      flush(rsv);


end;
procedure riteracii(biskl:boolean);
const chi=75;
var i,ih,ih2:integer;
    lp,rog2,ad,bd,h,rog0,tetag0:extended;
    rsx:text;
    ti0,tiro,titeta,tirot,titetat,tir,tivr,tim,tip,vgrad:das;
    tiro0,titeta0,tirot0,titetat0,tir0,tivr0,tim0,tip0:das;

    ch:char;
    f0,f1,f2,f3,f4,maxnev,ff1,ff2:extended;
    mi,itr:integer;
    sbvr,sbmassa,sbpia:boolean;
    spia,svmassa,svvr,raz,ln:extended;
    uge:boolean;
    siro0,siteta0,srot,stetat,sr,eps,eps2,svll,svlll,min,
    fd1,fd2,fd3,fd4,fd5,fd6,fd7,fd8,mnt,tkmjd:extended;
    rsh,rsi:text;
    ir,mir,vht,nsxf,tsxf,irn:integer;
    bsxf,bvht:boolean;
    sxf:array[1..8] of ^ttrs;
    ved,mnved:stb;
    svgrad:das;




begin
      bklpvd:=true;

      btvht:=false;
      nevyh:=false;
      for i:=1 to 8 do new(sxf[i]);

      assign(rsi,'allrsi.res');
      rewrite(rsi);
      {eps:=1e-10;
      eps2:=1e-9;
      eps:= eps /1000000;
      eps2:=eps2/1000000;}
      eps:=1e-7;
      eps2:=1e-6;


      {l:=le;}


      siro0:=iro0;
      siteta0:=iteta0;
      srot:=rot;
      stetat:=tetat;
      sr:=r;

      uci:=false;
      ln:=l;
      uge:=false;

      svmassa:=massa;
      svvr:=vr;

      sbvr:=bvr;
      sbmassa:=bmassa;
      sbpia:=bpia;
      spia:=pia;

      {bvr:=false;
      bmassa:=false;
      bpia:=false;
      pia:=1;}
      {br:=false;
      bro0:=false;
      bteta0:=false;}
      {brot:=false;
      btetat:=false;}


      dro:=eps*180/pi;
      ch:='p';
         vht:=0;
         fd:=0;
         bsxf:=false;
         nsxf:=0;

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
        if bvr then mpdq.tpvd(iro0,iteta0,rot,tetat,r,(vr-0.5*eps),massa,et^[1],pia,tivr0);
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
                          tcir2(ti0,tiro0,tiro,titeta0,titeta,
                               tirot,tirot0,titetat0,titetat,
                                tir0,tir,tivr0,tivr,tim0,tim,tip0,tip,
                                rog0,tetag0,eps,raz,i,
                                jb,id,
                                et,ro,teta,
                                rog,tetag,rog1,tetag1,dro,
                                bro0,bteta0,brot,btetat,br,bvr,bmassa,bpia,
                                al
                                 );
                     end;
                end;
                {halt(0);}
                {close(rsx);}
                ptlst4(al,xf,cv,i*2,jb);
                nsxf:=nsxf+1;
                if nsxf=(trs+1) then
                begin
                     nsxf:=1;
                     bsxf:=true;
                end;
                for ih:=1 to jb do
                begin
                     sxf[ih]^[nsxf]:=abs(xf[ih]);
                end;
                {sgmnk(a,cv,xf,sgm,jb,i*2);}
                {write(rsi,fd:5);
                for ih:=1 to jb do
                begin
                     write(rsi,' ',abs(xf[ih]):20:15);
                end;
                writeln(rsi);}

                for ih:=1 to jb do vgrad[ih]:=xf[ih];


                f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
                if f0>1e37 then
                begin
                     gf0:=f0;
                end;
                writeln(rsi,fd:5,' ',f0:20:15,' ',vht:4);
                gf0:=f0;
                gprd:=ti0[4];
                gnak:=ti0[11];



                f2:=flam(l,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                l:=abs(l);
                svll:=l;
                while (f2>f0) and (l>1e-10) do
                begin
                     l:=l/1.1;
                     f2:=flam(l,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                end;
                bvht:=(l<1e-10);
                if (l<1e-10) {and false} {and( svll>1e-7)} {or (abs(ev-1)<0.01)} then
                begin


                     tsxf:=nsxf;
                     if bsxf then tsxf:=trs;
                     for ih:=1 to jb do
                     begin
                          sort3(sxf[ih]^,tsxf);
                          ih2:=trunc(tsxf/2);
                          if ih2<1 then ih2:=1;
                          if ih2>trs then ih2:=trs;


                          svgrad[ih]:=sxf[ih]^[ih2];
                     end;
                     {svgrad[2]:=0;}
                     min:=1e+38;
                     lp:=1e-5;
                     if svll>lp then lp:=svll;
                     {irn:=0;
                     repeat}
                           ir:=0;
                           repeat
                           begin
                                ir:=ir+1;
                                for ih:=1 to jb do ved[ih]:=2*random-1;
                                for ih:=1 to jb do
                                begin
                                     vgrad[ih]:=svgrad[ih]*ved[ih];
                                     {vgrad[ih]:=xf[ih]*ved[ih];}
                                end;
                                fd8:=flam(lp,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;
                                if fd8<min then
                                begin
                                     min:=fd8;
                                     for ih:=1 to jb do mnved[ih]:=ved[ih];
                                end;
                           end;
                           until ((min<0) and (ir>chi)) or (ir>chi*10);
                           {if min>=0 then
                           begin
                                lp:=lp/10;
                           end;}
                           {irn:=irn+1;
                     until (min<0) or (irn>10);}
                     if min<0 then
                     begin
                          mnt:=1;
                          for ih:=1 to jb do
                          begin
                               xf[ih]:=svgrad[ih]*mnved[ih];
                               {xf[ih]:=xf[ih]*mnved[ih];}

                          end;

                          for ih:=1 to jb do vgrad[ih]:=xf[ih];
                          fd8:=min;

                          repeat
                                fd7:=fd8;
                                mnt:=mnt*2;
                                fd8:=flam(lp*mnt,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;
                          until (fd7<fd8) or (mnt*lp>1);
                          mnt:=mnt/2;
                          l:=lp*mnt;
                          {l:=lp;}
                          vht:=vht+1;
                     end;

                     fd1:=0;
                     {for ih:=1 to jb do vgrad[ih]:=xf[ih];
                     vgrad[1]:=0;
                     fd1:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;

                     for ih:=1 to jb do vgrad[ih]:=xf[ih];
                     vgrad[2]:=0;
                     fd2:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;

                     for ih:=1 to jb do vgrad[ih]:=xf[ih];
                     vgrad[3]:=0;
                     fd3:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;

                     for ih:=1 to jb do vgrad[ih]:=xf[ih];
                     vgrad[4]:=0;
                     fd4:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;

                     for ih:=1 to jb do vgrad[ih]:=xf[ih];
                     vgrad[5]:=0;
                     fd5:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;

                     for ih:=1 to jb do vgrad[ih]:=xf[ih];
                     vgrad[6]:=0;
                     fd6:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;

                     for ih:=1 to jb do vgrad[ih]:=xf[ih];
                     vgrad[7]:=0;
                     fd7:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;

                     for ih:=1 to jb do vgrad[ih]:=xf[ih];
                     vgrad[1]:=0;
                     vgrad[5]:=0;
                     vgrad[7]:=0;
                     fd8:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;
                     f0:=0;}

                     (*min:=1e+38;
                     lp:=1e-5;
                     if svll>lp then lp:=svll;
                     for ir:=1 to jb do
                     begin
                          for ih:=1 to jb do vgrad[ih]:=xf[ih];
                          vgrad[ir]:=0;
                          fd8:=flam(lp,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;
                          if fd8<min then
                          begin
                               min:=fd8;
                               mir:=ir;
                          end;
                     end;
                     if min<0 then
                     begin
                          xf[mir]:=0;
                          l:=lp;
                          mnt:=1;
                          for ih:=1 to jb do vgrad[ih]:=xf[ih];

                          repeat
                                fd7:=fd8;
                                mnt:=mnt*2;
                                fd8:=flam(lp*mnt,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;
                          until (fd7<fd8) or (mnt*lp>1);
                          mnt:=mnt/2;
                          l:=lp*mnt;
                          vht:=vht+1;
                     end else
                     begin

                          min:=1e+38;
                          lp:=1e-5/t1024;
                          if svll>lp then lp:=svll;
                          for ir:=1 to jb do
                          begin
                              for ih:=1 to jb do vgrad[ih]:=0;
                              vgrad[ir]:=xf[ir];
                              fd8:=flam(lp,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;
                              if fd8<min then
                              begin
                                   min:=fd8;
                                   mir:=ir;
                              end;
                          end;
                          if min<0 then
                          begin
                              for ih:=1 to jb do
                              begin
                                  if ih<>mir then xf[ih]:=0;
                              end;
                              for ih:=1 to jb do vgrad[ih]:=xf[ih];
                               l:=lp;
                               mnt:=1;
                               repeat
                                     fd7:=fd8;
                                     mnt:=mnt*2;
                                     fd8:=flam(lp*mnt,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;
                               until (fd7<fd8) or (mnt*lp>1);
                               mnt:=mnt/2;
                               l:=lp*mnt;
                               vht:=vht+1;
                           end else
                           begin
                                {svll:=-1e-5;
                                for ih:=1 to jb do vgrad[ih]:=0;
                                vgrad[1]:=xf[1];
                                fd1:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;

                                for ih:=1 to jb do vgrad[ih]:=0;
                                vgrad[2]:=xf[2];
                                fd2:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;

                                for ih:=1 to jb do vgrad[ih]:=0;
                                vgrad[3]:=xf[3];
                                fd3:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;

                                for ih:=1 to jb do vgrad[ih]:=0;
                                vgrad[4]:=xf[4];
                                fd4:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;

                                for ih:=1 to jb do vgrad[ih]:=0;
                                vgrad[5]:=xf[5];
                                fd5:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;

                                for ih:=1 to jb do vgrad[ih]:=0;
                                vgrad[6]:=xf[6];
                                fd6:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;

                                for ih:=1 to jb do vgrad[ih]:=0;
                                vgrad[7]:=xf[7];
                                fd7:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;

                                for ih:=1 to jb do vgrad[ih]:=0;
                                vgrad[8]:=xf[8];
                                fd8:=flam(svll,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;
                                f0:=0;}

                                min:=1e+38;
                                lp:=1e-5/t1024;
                                if svll>lp then lp:=svll;
                                for ir:=1 to jb do
                                begin
                                     for ih:=1 to jb do vgrad[ih]:=0;
                                     vgrad[ir]:=-xf[ir];
                                     fd8:=flam(lp,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;
                                     if fd8<min then
                                     begin
                                          min:=fd8;
                                          mir:=ir;
                                     end;
                                end;
                                if min<0 then
                                begin
                                     for ih:=1 to jb do
                                     begin
                                         if ih<>mir then xf[ih]:=0 else xf[ih]:=-xf[ih];
                                     end;
                                     for ih:=1 to jb do vgrad[ih]:=xf[ih];
                                     l:=lp;
                                     mnt:=1;
                                     repeat
                                           fd7:=fd8;
                                           mnt:=mnt*2;
                                           fd8:=flam(lp*mnt,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;
                                     until (fd7<fd8) or (mnt*lp>1);
                                     mnt:=mnt/2;
                                     l:=lp*mnt;
                                     vht:=vht+1;
                                end;

                           end;




                     end;*)

                end;
                (*if (l<1e-10) {and( svll>1e-7)} then
                begin
                     svlll:=l;
                     l:=svll;
                     f2:=flam(l,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                     for ih:=1 to jb do vgrad[ih]:=-xf[ih];
                     while (f2>f0) {and (fd<>258)} do
                     begin
                          l:=l/1.1;
                          f2:=flam(l,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                     end;
                     if l>1e-10 then
                     begin
                          for ih:=1 to jb do xf[ih]:=-xf[ih];
                     end else
                     begin
                          l:=svlll;
                     end;

                end;*)


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
                if not(bvht) then
                begin
                     f3:=flam(l/1.1,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                     f4:=flam(l*1.03,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);
                     if (f3<f2) and (f3<f0) then l:=l/1.1;
                     if (f4<=f2) and (f4<f0) then l:=l*1.03;
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
                if massa<0 then
                begin
                     bvr:=sbvr;
                     bmassa:=sbmassa;
                     bpia:=sbpia;
                     pia:=spia;
                     gf0:=1e38;

                     exit;

                end;



                if (fd mod 1 =0) and not(bvnevyd) then
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


                writeln('iteracija=',fd:5,' kol-vo urav=',i:5,' vht=',vht:2 );

                writeln('f0(oev)=',f0:10:5 );
                writeln('lp=',lp );
                writeln('e=',ev:10:5 );
                tkmjd:=tekmjdr-bmjd;
                tkmjd:=tkmjd*60*24;
                writeln('tkmjd=',tkmjd:10:5 );


                {readln;}
                if not(uge) then
                begin
                     ku2:=ku2+1;
                     uge:=true;
                end;

                if biskl then wvh;
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

            {if (l<1e-10) or (ev>2) or (abs(tetat)>20) then
            begin
                 vyh:=true;
            end;}
            if fd=365 then
            begin
                  ch:='r';
            end;
            if ch='a' then halt(0);
            tkmjd:=tekmjdr-bmjd;
            tkmjd:=tkmjd*60*24;


            nevyh:={(abs(ev-1.0)<0.001) or}
                   {(ev>2) or}
                   (r>3000) or
                   (fd>t50000) or
                   (abs(tetat)>100) or
                   {((fd>200) and ((vht*10)>fd)) or
                   ((fd>10) and ((vht*2)>fd)) or}
                   (massa>1000) or
                   (massa<0) or
                   ((abs(pia-1)<0.001) and (abs(vr)>10)) or
                   (abs(rot)>10) or
                   (l<1e-10) or
                   (tkmjd>predmin) or
                   (ti0[4]<5)
                    ;
            if koi then
            begin
                   nevyh:=nevyh or ((fd>200) and ((vht*10)>fd)) or
                   ((fd>10) and ((vht*2)>fd));


            end;
            if nevyh and (fd>200) then
            begin
                  ch:='r';
            end;



          until (ch='a') or vyh or nevyh;
          if {(l>1e-10) and} (ev<50) then
          begin

          if biskl then
          begin
               if (abs(ev-1)<0.001) and (t50000<10000) then
               begin
                    fd:=0;
               end;
               maxiskl(iro0,iteta0,rot,tetat,r,vr,massa,pia,oevr,oevt,maxnev,mi);
               ew:=false;
               begin
                    ew:=maxnev>porog;
                    if ew then
                    begin
                         cvb2(mi,3);
                         {if nevyh then
                         begin
                               iro0:=siro0;
                               iteta0:=siteta0;
                               rot:=srot;
                               tetat:=stetat;
                               r:=sr;
                               massa:=svmassa;
                               vr:=svvr;

                         end;}
                         l:=ln;

                    end;
               end;
          end
          else ew:=false;
          end;
          wvh;


          until  (ch='a') or not(ew) or (tkmjd>predmin) ;

          {if biskl then vydres2(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt,oevr,oevt);}
          end;
          uci:=true;
          if ch='a' then uci:=false;
          {if l<1e-10 then uci:=false;}
          if ev>2 then uci:=false;
          if (r>3000) or (fd>t50000) or (abs(tetat)>1000) then uci:=false;
          if nevyh then uci:=false;

          uci:=true;

          if uci then
          begin
                sgmnk2(al,cv,xf,sgm,jb,i*2);

          end else
          begin
               gf0:=1e39;
                               iro0:=siro0;
                               iteta0:=siteta0;
                               rot:=srot;
                               tetat:=stetat;
                               r:=sr;
                               massa:=svmassa;
                               vr:=svvr;


          end;
      {massa:=svmassa;
      vr:=svvr;}

      bvr:=sbvr;
      bmassa:=sbmassa;
      bpia:=sbpia;
      pia:=spia;
      close(rsi);
      for i:=1 to 8 do dispose(sxf[i]);
      btvht:=((fd>200) and ((vht*10)>fd)) or ((fd>10) and ((vht*2)>fd));




end;

procedure getprelem2(ev,th,prd:extended;var ti:das);
var aa,bb,ff,gg,oev:extended;
begin
                 GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);

                 ti[1]:=0;
                 ti[2]:=ev;

                 ti[3]:=th;
                 ti[4]:=prd;

                 ti[5]:=aa;
                 ti[6]:=bb;
                 ti[7]:=ff;
                 ti[8]:=gg;

                 ti[9]:=0;
                 ti[10]:=0;
                 ti[11]:=0;
                 ti[12]:=0;

end;

procedure  getelem(abp,ev,th,prd,ai,omm,omb,t,pia:extended;var ti:das);
var comm,somm,comb,somb,ci,si,aa,bb,ff,gg:extended;
begin
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


                 ti[1]:=abp;
                 ti[2]:=ev;

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

end;
procedure tcir3(
               {ti0,tiro0,tiro,titeta0,titeta,
               tirot,tirot0,titetat0,titetat,
               tir0,tir,tivr0,tivr,tim0,tim
               :das;}
               ti0,tiabp0,tiabp,tiev0,tiev,
               tith0,tith,tiprd0,tiprd,
               tiai0,tiai,tiomm0,tiomm,tiomb0,tiomb
               :das;
               var rog0,tetag0,eps,raz:extended;
               var i:integer;
               var jb,id:integer;
               et,ro,teta:dst3;
               var rog,tetag,rog1,tetag1,dro:extended;

               a:msls2);
var
    bro0,bteta0,brot,btetat,br,bvr,bmassa,bpia:extended;
    iro0,iteta0,rot,tetat,r,vr,massa:boolean;

begin
                     jb:=0;
                     {i:=id;}
                     i:=i+1;
                     mpdq.rtti(ti0,et^[id],rog,tetag);
                     {abp,ev,th,prd,ai,omm,omb}
                     if uabp then
                     begin
                          mpdq.rtti(tiabp0,et^[id],rog0,tetag0);
                          mpdq.rtti(tiabp,et^[id],rog1,tetag1);
                          a^[i*2-1]^[jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;

                     if uev then
                     begin
                          mpdq.rtti(tiev0,et^[id],rog0,tetag0);
                          mpdq.rtti(tiev,et^[id],rog1,tetag1);

                          a^[i*2-1]^[jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;
                     if uth then
                     begin
                          mpdq.rtti(tith0,et^[id],rog0,tetag0);
                          mpdq.rtti(tith,et^[id],rog1,tetag1);
                          a^[i*2-1]^[jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;
                     if uprd then
                     begin
                          mpdq.rtti(tiprd0,et^[id],rog0,tetag0);
                          mpdq.rtti(tiprd,et^[id],rog1,tetag1);
                          a^[i*2-1]^[jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;
                     if uai then
                     begin
                          mpdq.rtti(tiai0,et^[id],rog0,tetag0);
                          mpdq.rtti(tiai,et^[id],rog1,tetag1);
                          a^[i*2-1]^[jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;


                     if uomm then
                     begin
                          mpdq.rtti(tiomm0,et^[id],rog0,tetag0);
                          mpdq.rtti(tiomm,et^[id],rog1,tetag1);
                          a^[i*2-1]^[jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;
                     if uomb then
                     begin
                          mpdq.rtti(tiomb0,et^[id],rog0,tetag0);
                          mpdq.rtti(tiomb,et^[id],rog1,tetag1);
                          a^[i*2-1]^[jb+1]:=(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;

                     a^[i*2-1]^[jb+1]:=rog-ro^[id];
                     raz:=(tetag-teta^[id]);
                     ogr2(raz);
                      a^[i*2]^[jb+1]:=raz*pi*ro^[id]/180;
end;

procedure tcir4(
               {ti0,tiro0,tiro,titeta0,titeta,
               tirot,tirot0,titetat0,titetat,
               tir0,tir,tivr0,tivr,tim0,tim
               :das;}
               ti0,tiabp0,tiabp,tiev0,tiev,
               tith0,tith,tiprd0,tiprd,
               tiai0,tiai,tiomm0,tiomm,tiomb0,tiomb
               :das;
               var rog0,tetag0,eps,raz:extended;
               var i:integer;
               var jb,id:integer;
               et,ro,teta:dst3;
               var rog,tetag,rog1,tetag1,dro:extended;

               a:msls2);
var
    bro0,bteta0,brot,btetat,br,bvr,bmassa,bpia:extended;
    iro0,iteta0,rot,tetat,r,vr,massa:boolean;

begin
                     jb:=0;
                     {i:=id;}
                     i:=i+1;
                     mpdq.rtti4(ti0,et^[id],rog,tetag);
                     {abp,ev,th,prd,ai,omm,omb}
                     if uev then
                     begin
                          mpdq.rtti4(tiev0,et^[id],rog0,tetag0);
                          mpdq.rtti4(tiev,et^[id],rog1,tetag1);

                          a^[i*2-1]^[jb+1]:=aves^[id]*(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=aves^[id]*raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;
                     if uth then
                     begin
                          mpdq.rtti4(tith0,et^[id],rog0,tetag0);
                          mpdq.rtti4(tith,et^[id],rog1,tetag1);
                          a^[i*2-1]^[jb+1]:=aves^[id]*(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=aves^[id]*raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;
                     if uprd then
                     begin
                          mpdq.rtti4(tiprd0,et^[id],rog0,tetag0);
                          mpdq.rtti4(tiprd,et^[id],rog1,tetag1);
                          a^[i*2-1]^[jb+1]:=aves^[id]*(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=aves^[id]*raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;
                     a^[i*2-1]^[jb+1]:=aves^[id]*(rog-ro^[id]);
                     raz:=(tetag-teta^[id]);
                     ogr2(raz);
                     a^[i*2]^[jb+1]:=aves^[id]*raz*pi*ro^[id]/180;
end;

function fnev3(abp,ev,th,prd,ai,omm,omb:extended):extended;
var
ti0:das;
dro,dte,raz,max:extended;
nt2:integer;
begin
      nevgmi:=-1;
      max:=-1e38;
      getelem(abp,ev,th,prd,ai,omm,omb,et^[1],pia,ti0);


           if (abs(ti0[2]-1.0)<0.000000001) or (ti0[2]>10) or (ti0[4]<5) or (ti0[2]<0) or (ev>0.99) then
           begin
                fnev3:=1e39;
                exit;

           end;

      sum:=0;
      nt2:=0;
      for id:=1 to nt do
      begin
           if gbop^[id]=0 then
           begin
                nt2:=nt2+1;
                mpdq.rtti(ti0,et^[id],rog,tetag);
                dro:=rog-ro^[id];
                raz:=(tetag-teta^[id]);
                ogr2(raz);
                dte:=raz*pi*ro^[id]/180;
                {dte:=raz*pi*rog/180;}
                sum:=sum+dro*dro+dte*dte;
                dro:=abs(dro);
                dte:=abs(dte);
                if dro>max then
                begin
                     max:=dro;
                     nevgmi:=id;
                end;
                if dte>max then
                begin
                     max:=dte;
                     nevgmi:=id;

                end;
           end;
      end;
      sum:=0.5*sum/nt2;
      sum:=sqrt(sum);
      fnev3:= sum;

end;

function fnev4(ev,th,prd:extended):extended;
var

dro,dte,raz,max:extended;
nt2:integer;
aa,bb,ff,gg,oev:extended;
begin
      nevgmi:=-1;
      max:=-1e38;
      GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
      if not((swds='14575-2125') and (comp10='Ba,Bb')) then
      begin
           if (abs(ev-1.0)<0.000000001) or (ev>10) or (prd<1) or (ev<0) or (ev>0.99) then
           begin
                fnev4:=1e38;
                exit;
           end;
      end else
      begin
           if (abs(ev-1.0)<0.000000001) or (ev>10) or (prd<0.01) or (ev<0) or (ev>0.99) then
           begin
                fnev4:=1e38;
                exit;
           end;

      end;

      sum:=0;
      nt2:=0;
      for id:=1 to nt do
      begin
           if gbop^[id]=0 then
           begin
                nt2:=nt2+1;
                mpdq.rtti3(et^[id],ev,th,prd,aa,bb,ff,gg,rog,tetag);
                dro:=rog-ro^[id];
                raz:=(tetag-teta^[id]);
                ogr2(raz);
                dte:=raz*pi*ro^[id]/180;
                {dte:=raz*pi*rog/180;}
                dro:=abs(dro);
                dte:=abs(dte);
                if dro>max then
                begin
                     max:=dro;
                     nevgmi:=id;
                end;
                if dte>max then
                begin
                     max:=dte;
                     nevgmi:=id;

                end;

                dro:=dro*aves^[id];
                dte:=dte*aves^[id];
                sum:=sum+dro*dro+dte*dte;
           end;
      end;
      sum:=0.5*sum/nt2;
      sum:=sqrt(sum);
      fnev4:= sum;

end;
function flam3(lam:extended;vgrad:das;abp,ev,th,prd,ai,omm,omb:extended):extended;
begin
                     jb:=0;
                     if uabp then
                     begin
                          jb:=jb+1;
                          abp:=abp-lam*vgrad[jb];
                     end;

                     if uev then
                     begin
                          jb:=jb+1;
                          ev:=ev-lam*vgrad[jb];
                     end;
                     if uth then
                     begin
                          jb:=jb+1;
                          th:=th-lam*vgrad[jb];

                     end;
                     {prd,ai,omm,omb}
                     if uprd then
                     begin
                          jb:=jb+1;
                          prd:=prd-lam*vgrad[jb];
                     end;
                     if uai then
                     begin
                          jb:=jb+1;
                          ai:=ai-lam*vgrad[jb];
                     end;


                     if uomm then
                     begin
                          jb:=jb+1;
                          omm:=omm-lam*vgrad[jb];
                     end;
                     if uomb then
                     begin
                          jb:=jb+1;
                          omb:=omb-lam*vgrad[jb];
                     end;


        flam3:=fnev3(abp,ev,th,prd,ai,omm,omb);

end;
function flam4(lam:extended;vgrad:das;ev,th,prd:extended):extended;
begin
                     jb:=0;
                     if uev then
                     begin
                          jb:=jb+1;
                          ev:=ev-lam*vgrad[jb];
                     end;
                     if uth then
                     begin
                          jb:=jb+1;
                          th:=th-lam*vgrad[jb];

                     end;
                     {prd,ai,omm,omb}
                     if uprd then
                     begin
                          jb:=jb+1;
                          prd:=prd-lam*vgrad[jb];
                     end;
        flam4:=fnev4(ev,th,prd);

end;

procedure mappend(var rs:text);
begin
{$I-}
     append(rs);
     if ioresult<>0 then rewrite(rs);
{$I+}
end;



procedure elemiter(biskl:boolean);
const chi=75;
var i,ih,ih2:integer;
    lp,rog2,ad,bd,h,rog0,tetag0,izmp:extended;
    rsx:text;
    ti0{,tiro,titeta,tirot,titetat,tir,tivr,tim,tip},vgrad:das;
    {tiro0,titeta0,tirot0,titetat0,tir0,tivr0,tim0,tip0:das;}
    tiabp,tiev,tith,tiprd,tiai,tiomm,tiomb:das;
    tiabp0,tiev0,tith0,tiprd0,tiai0,tiomm0,tiomb0:das;

    ch:char;
    f0,f1,f2,f3,f4,maxnev,ff1,ff2:extended;
    mi,itr:integer;
    sbvr,sbmassa,sbpia:boolean;
    spia,svmassa,svvr,raz,ln:extended;
    uge:boolean;
    siro0,siteta0,srot,stetat,sr,eps,eps2,svll,svlll,min,
    fd1,fd2,fd3,fd4,fd5,fd6,fd7,fd8,mnt,tkmjd:extended;
    rsh,rsi:text;
    ir,mir,vht,nsxf,tsxf,irn:integer;
    bsxf,bvht:boolean;
    sxf:array[1..8] of ^ttrs;
    ved,mnved:stb;
    svgrad:das;
    rso2:text;
    {bro0,bteta0,brot,btetat,br,bvr,bmassa,bpia:extended;}
    {iro0,iteta0,rot,tetat,r,vr,massa:boolean;}






begin
     bklpvd:=false;

     {assign(rso2,'kaiomb.res');
     mappend(rso2);
     writeln(rso2,ai*180/pi:20:15,' ',omb*180/pi:20:15);
     close(rso2);}
     ClrScr;
     {l:=1e-5;}
     {mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);}
     {mpdq.tpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);
     abp:=ti0[1];
     ev:= ti0[2];
     th:= ti0[3];
     prd:=ti0[4];
     omm:=ti0[9];
     omb:=ti0[10];
     ai:= ti0[11];}


     {f0:=fnev3(abp,ev,th,prd,ai,omm,omb);
     mbcir:=true;
     f2:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
     mbcir:=false;
     f3:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);}
     {uth:=false;}

     {ev:=0.995;
     uev:=false;}
     {uprd:=false;}

      btvht:=false;
      nevyh:=false;
      for i:=1 to 8 do new(sxf[i]);
      eps:=1e-7;
      eps2:=1e-6;
      uci:=false;
      ln:=l;
      uge:=false;
      dro:=eps*180/pi;
        ch:='p';
         vht:=0;
         fd:=0;
         bsxf:=false;
         nsxf:=0;

          begin
                {massa:=massa-0.1;}
          repeat
          t3:=nt;

          {uabp,uev,uth,uprd,uai,uomm,uomb}

          repeat
                i:=0;
                getelem(abp,ev,th,prd,ai,omm,omb,et^[1],pia,ti0);
                ev:=ti0[2];
        if uabp then getelem(abp-0.5*eps,ev,th,prd,ai,omm,omb,et^[1],pia,tiabp0);
        if uev then  getelem(abp,ev-0.5*eps,th,prd,ai,omm,omb,et^[1],pia,tiev0);
        if uth then  getelem(abp,ev,th-0.5*eps,prd,ai,omm,omb,et^[1],pia,tith0);
        if uprd then  getelem(abp,ev,th,prd-0.5*eps,ai,omm,omb,et^[1],pia,tiprd0);
        if uai then  getelem(abp,ev,th,prd,ai-0.5*eps,omm,omb,et^[1],pia,tiai0);
        if uomm then  getelem(abp,ev,th,prd,ai,omm-0.5*eps,omb,et^[1],pia,tiomm0);
        if uomb then  getelem(abp,ev,th,prd,ai,omm,omb-0.5*eps,et^[1],pia,tiomb0);

        if uabp then getelem(abp+0.5*eps,ev,th,prd,ai,omm,omb,et^[1],pia,tiabp);
        if uev then  getelem(abp,ev+0.5*eps,th,prd,ai,omm,omb,et^[1],pia,tiev);
        if uth then  getelem(abp,ev,th+0.5*eps,prd,ai,omm,omb,et^[1],pia,tith);
        if uprd then  getelem(abp,ev,th,prd+0.5*eps,ai,omm,omb,et^[1],pia,tiprd);
        if uai then  getelem(abp,ev,th,prd,ai+0.5*eps,omm,omb,et^[1],pia,tiai);
        if uomm then  getelem(abp,ev,th,prd,ai,omm+0.5*eps,omb,et^[1],pia,tiomm);
        if uomb then  getelem(abp,ev,th,prd,ai,omm,omb+0.5*eps,et^[1],pia,tiomb);

                i:=0;
                for id:=1 to t3 do
                begin
                     if gbop^[id]=0 then
                     begin
                          {tcir2(ti0,tiro0,tiro,titeta0,titeta,
                               tirot,tirot0,titetat0,titetat,
                                tir0,tir,tivr0,tivr,tim0,tim,tip0,tip,
                                rog0,tetag0,eps,raz,i,
                                jb,id,
                                et,ro,teta,
                                rog,tetag,rog1,tetag1,dro,
                                bro0,bteta0,brot,btetat,br,bvr,bmassa,bpia,
                                al
                                 );}
                          {abp,ev,th,prd,ai,omm,omb}
                          tcir3(ti0,tiabp0,tiabp,tiev0,tiev,
                               tith0,tith,tiprd0,tiprd,
                                tiai0,tiai,tiomm0,tiomm,tiomb0,tiomb,
                                rog0,tetag0,eps,raz,i,
                                jb,id,
                                et,ro,teta,
                                rog,tetag,rog1,tetag1,dro,
                                al
                                 );
                     end;
                end;
                ptlst4(al,xf,cv,i*2,jb);
                {for ih:=1 to jb do xf[ih]:=-xf[ih];}


                nsxf:=nsxf+1;
                if nsxf=(trs+1) then
                begin
                     nsxf:=1;
                     bsxf:=true;
                end;
                for ih:=1 to jb do
                begin
                     sxf[ih]^[nsxf]:=abs(xf[ih]);
                end;
                for ih:=1 to jb do vgrad[ih]:=xf[ih];


                f0:=fnev3(abp,ev,th,prd,ai,omm,omb);
                if f0>1e37 then
                begin
                     gf0:=f0;
                end;
                gf0:=f0;
                gprd:=ti0[4];
                gnak:=ti0[11];



                f2:=flam3(l,vgrad,abp,ev,th,prd,ai,omm,omb);
                l:=abs(l);
                svll:=l;
                while (f2>f0) and (l>1e-10) do
                begin
                     l:=l/10;
                     f2:=flam3(l,vgrad,abp,ev,th,prd,ai,omm,omb);
                end;
                bvht:=(l<1e-10);
                if (l<1e-10) {and false} {and( svll>1e-7)} or (abs(ev-1)<0.01) then
                begin


                     tsxf:=nsxf;
                     if bsxf then tsxf:=trs;
                     for ih:=1 to jb do
                     begin
                          sort3(sxf[ih]^,tsxf);
                          ih2:=trunc(tsxf/2);
                          if ih2<1 then ih2:=1;
                          if ih2>trs then ih2:=trs;


                          svgrad[ih]:=sxf[ih]^[ih2];
                     end;
                     {svgrad[2]:=0;}
                     min:=1e+38;
                     lp:=1e-5;
                     if svll>lp then lp:=svll;
                     {irn:=0;
                     repeat}
                           ir:=0;
                           repeat
                           begin
                                ir:=ir+1;
                                for ih:=1 to jb do ved[ih]:=2*random-1;
                                for ih:=1 to jb do
                                begin
                                     vgrad[ih]:=svgrad[ih]*ved[ih];
                                     {vgrad[ih]:=xf[ih]*ved[ih];}
                                end;
                                {fd8:=flam(lp,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;}
                                fd8:=flam3(lp,vgrad,abp,ev,th,prd,ai,omm,omb)-f0;
                                if fd8<min then
                                begin
                                     min:=fd8;
                                     for ih:=1 to jb do mnved[ih]:=ved[ih];
                                end;
                           end;
                           until ((min<0) and (ir>chi)) or (ir>chi*3);
                           {if min>=0 then
                           begin
                                lp:=lp/10;
                           end;}
                           {irn:=irn+1;
                     until (min<0) or (irn>10);}
                     if min<0 then
                     begin
                          mnt:=1;
                          for ih:=1 to jb do
                          begin
                               xf[ih]:=svgrad[ih]*mnved[ih];
                               {xf[ih]:=xf[ih]*mnved[ih];}

                          end;

                          for ih:=1 to jb do vgrad[ih]:=xf[ih];
                          fd8:=min;

                          repeat
                                fd7:=fd8;
                                mnt:=mnt*2;
                                {fd8:=flam(lp*mnt,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;}
                                fd8:=flam3(lp*mnt,vgrad,abp,ev,th,prd,ai,omm,omb)-f0;
                          until (fd7<fd8) or (mnt*lp>1);
                          mnt:=mnt/2;
                          l:=lp*mnt;
                          vht:=vht+1;
                     end;

                     fd1:=0;
                end;
                if not(bvht) then
                begin
                     {f3:=flam(l/1.1,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);}
                     f3:=flam3(l/1.1,vgrad,abp,ev,th,prd,ai,omm,omb);
                     {f4:=flam(l*1.03,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);}
                     f4:=flam3(l*1.03,vgrad,abp,ev,th,prd,ai,omm,omb);
                     if (f3<f2) and (f3<f0) then l:=l/1.1;
                     if (f4<=f2) and (f4<f0) then l:=l*1.03;
                end;
                {l:=1e-2;}
                {if fd>1000 then l:=1e-1;}


                     jb:=0;
                     vyh:=true;

                     lp:=l;
                     if uabp then
                     begin
                          jb:=jb+1;
                          abp:=abp-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<eps2);
                     end;

                     if uev then
                     begin
                          jb:=jb+1;
                          ev:=ev-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<eps2);
                     end;
                     if uth then
                     begin

                          jb:=jb+1;
                          th:=th-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<eps2);
                     end;
                     {prd,ai,omm,omb}
                     if uprd then
                     begin
                          jb:=jb+1;
                          prd:=prd-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<eps2);
                          izmp:=-lp*xf[jb];

                     end;
                     if uai then
                     begin
                          jb:=jb+1;
                          ai:=ai-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<eps2);

                     end;

                     if uomm then
                     begin
                          jb:=jb+1;;
                          omm:=omm-lp*xf[jb];

                          vyh:=vyh and (abs(xf[jb])<eps2);

                     end;

                     if uomb then
                     begin
                          jb:=jb+1;
                          omb:=omb-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<eps2);

                     end;

                fd:=fd+1;
                {if massa<0 then
                begin
                     bvr:=sbvr;
                     bmassa:=sbmassa;
                     bpia:=sbpia;
                     pia:=spia;
                     gf0:=1e39;

                     exit;

                end;}



                if (fd mod 1 =0) and not(bvnevyd) then
                begin
                {abp,ev,th,prd,ai,omm,omb}

                gotoxy(1,1);
                writeln;
                if uabp then pst:='     var'
                        else pst:=' not var'  ;
                writeln('   abp=',abp:10:5,pst);
                if uev  then pst:='     var'
                        else pst:=' not var'  ;
                writeln('    ev=',ev:10:5,pst);
                if uth  then pst:='     var'
                        else pst:=' not var'  ;
                writeln('    th=',th+t0 :10:5,pst);
                if uprd then pst:='     var'
                        else pst:=' not var'  ;
                writeln('   prd=',prd  :10:5,pst);
                if uai then pst:='     var'
                       else pst:=' not var'  ;
                writeln('     ai=',ai*180/pi:10:5,pst);
                if uomm then pst:='     var'
                        else pst:=' not var'  ;
                writeln('    omm=',omm*180/pi:10:5,pst);
                if uomb   then pst:='     var'
                          else pst:=' not var'  ;
                writeln('    omb=',omb*180/pi:10:5,pst);


                writeln('iteracija=',fd:5,' kol-vo urav=',i:5,' vht=',vht:2 );

                writeln('f0(oev)=',f0:15:20 );
                writeln('lp=',lp );
                writeln('e=',ev:10:5 );
                tkmjd:=tekmjdr-bmjd;
                tkmjd:=tkmjd*60*24;
                writeln('tkmjd=',tkmjd:10:5 );
                end;


                if not(uge) then
                begin
                     ku2:=ku2+1;
                     uge:=true;
                end;

                {if biskl then wvh;}
             ch:='p';

            if keypressed then ch:=readkey;
            if ch='a' then halt(0);
            tkmjd:=tekmjdr-bmjd;
            tkmjd:=tkmjd*60*24;


            nevyh:={(abs(ev-1.0)<0.001) or}
                   {(ev>2) or}
                   {(r>3000) or}
                   (fd>t50000) or
                   {(abs(tetat)>100) or}
                   {((fd>200) and ((vht*10)>fd)) or
                   ((fd>10) and ((vht*2)>fd)) or}
                   {(massa>1000) or
                   (massa<0) or
                   ((abs(pia-1)<0.001) and (abs(vr)>10)) or
                   (abs(rot)>10) or}
                   (l<1e-10) or
                   (tkmjd>predmin) or
                   (ti0[4]<5) {or}
                   {((prd>100000) and (izmp>0))}
                    ;
            if koi then
            begin
                   {nevyh:=nevyh or ((fd>200) and ((vht*10)>fd)) or
                   ((fd>10) and ((vht*2)>fd));}


            end;



          until (ch='a') or vyh or nevyh;

          if {(l>1e-10) and} (ev<50) then
          begin

          if biskl then
          begin
               maxiskl(iro0,iteta0,rot,tetat,r,vr,massa,pia,oevr,oevt,maxnev,mi);
               ew:=false;
               begin
                    ew:=maxnev>porog;
                    if ew then
                    begin
                         cvb2(mi,3);
                         {if nevyh then
                         begin
                               iro0:=siro0;
                               iteta0:=siteta0;
                               rot:=srot;
                               tetat:=stetat;
                               r:=sr;
                               massa:=svmassa;
                               vr:=svvr;

                         end;}
                         l:=ln;

                    end;
               end;
          end
          else ew:=false;
          end;
          wvh;


          until  (ch='a') or not(ew) or (tkmjd>predmin) ;

          {if biskl then vydres2(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt,oevr,oevt);}
          end;
          uci:=true;
          if ch='a' then uci:=false;
          {if l<1e-10 then uci:=false;}
          if ev>2 then uci:=false;
          if (r>3000) or (fd>t50000) or (abs(tetat)>1000) then uci:=false;
          if nevyh then uci:=false;

          uci:=true;

          if uci then
          begin
                sgmnk2(al,cv,xf,sgm,jb,i*2);

          end else
          begin
               gf0:=1e39;
                               iro0:=siro0;
                               iteta0:=siteta0;
                               rot:=srot;
                               tetat:=stetat;
                               r:=sr;
                               massa:=svmassa;
                               vr:=svvr;


          end;
      {massa:=svmassa;
      vr:=svvr;}

      for i:=1 to 8 do dispose(sxf[i]);
      btvht:=((fd>200) and ((vht*10)>fd)) or ((fd>10) and ((vht*2)>fd));




end;

procedure elemiter2(biskl:boolean);
const chi=75;
var i,ih,ih2:integer;
    lp,rog2,ad,bd,h,rog0,tetag0:extended;
    rsx:text;
    ti0{,tiro,titeta,tirot,titetat,tir,tivr,tim,tip},vgrad:das;
    {tiro0,titeta0,tirot0,titetat0,tir0,tivr0,tim0,tip0:das;}
    tiabp,tiev,tith,tiprd,tiai,tiomm,tiomb:das;
    tiabp0,tiev0,tith0,tiprd0,tiai0,tiomm0,tiomb0:das;

    ch:char;
    f0,f1,f2,f3,f4,maxnev,ff1,ff2:extended;
    mi,itr:integer;
    sbvr,sbmassa,sbpia:boolean;
    spia,svmassa,svvr,raz,ln:extended;
    uge:boolean;
    siro0,siteta0,srot,stetat,sr,eps,eps2,svll,svlll,min,
    fd1,fd2,fd3,fd4,fd5,fd6,fd7,fd8,mnt,tkmjd:extended;
    rsh,rsi:text;
    ir,mir,vht,nsxf,tsxf,irn:integer;
    bsxf,bvht:boolean;
    sxf:array[1..8] of ^ttrs;
    ved,mnved:stb;
    svgrad:das;
    {bro0,bteta0,brot,btetat,br,bvr,bmassa,bpia:extended;}
    {iro0,iteta0,rot,tetat,r,vr,massa:boolean;}
    qHour, qMinute, qSecond, qSec100,qqsecond: Word;






begin
     bklpvd:=false;

     ClrScr;
     l:=1e-5;
     {EV:=0.2;}


     f0:=fnev4(ev,th,prd);
     {uth:=false;}

     {ev:=0.995;
     uev:=false;}
     {uprd:=false;}

      btvht:=false;
      nevyh:=false;
      for i:=1 to 8 do new(sxf[i]);
      eps:=1e-7;
      eps2:=1e-6;
      uci:=false;
      ln:=l;
      uge:=false;
      dro:=eps*180/pi;
        ch:='p';
         vht:=0;
         fd:=0;
         bsxf:=false;
         nsxf:=0;

          begin
                {massa:=massa-0.1;}
          repeat
          t3:=nt;

          {uabp,uev,uth,uprd,uai,uomm,uomb}

          repeat
                i:=0;
                getprelem2(ev,th,prd,ti0);
                ev:=ti0[2];
        if uev then  getprelem2(ev-0.5*eps,th,prd,tiev0);
        if uth then  getprelem2(ev,th-0.5*eps,prd,tith0);
        if uprd then  getprelem2(ev,th,prd-0.5*eps,tiprd0);

        if uev then  getprelem2(ev+0.5*eps,th,prd,tiev);
        if uth then  getprelem2(ev,th+0.5*eps,prd,tith);
        if uprd then  getprelem2(ev,th,prd+0.5*eps,tiprd);
                i:=0;
                for id:=1 to t3 do
                begin
                     if gbop^[id]=0 then
                     begin
                          {tcir2(ti0,tiro0,tiro,titeta0,titeta,
                               tirot,tirot0,titetat0,titetat,
                                tir0,tir,tivr0,tivr,tim0,tim,tip0,tip,
                                rog0,tetag0,eps,raz,i,
                                jb,id,
                                et,ro,teta,
                                rog,tetag,rog1,tetag1,dro,
                                bro0,bteta0,brot,btetat,br,bvr,bmassa,bpia,
                                al
                                 );}
                          {abp,ev,th,prd,ai,omm,omb}
                          tcir4(ti0,tiabp0,tiabp,tiev0,tiev,
                               tith0,tith,tiprd0,tiprd,
                                tiai0,tiai,tiomm0,tiomm,tiomb0,tiomb,
                                rog0,tetag0,eps,raz,i,
                                jb,id,
                                et,ro,teta,
                                rog,tetag,rog1,tetag1,dro,
                                al
                                 );
                     end;
                end;
                ptlst4(al,xf,cv,i*2,jb);
                {for ih:=1 to jb do xf[ih]:=-xf[ih];}


                nsxf:=nsxf+1;
                if nsxf=(trs+1) then
                begin
                     nsxf:=1;
                     bsxf:=true;
                end;
                for ih:=1 to jb do
                begin
                     sxf[ih]^[nsxf]:=abs(xf[ih]);
                end;
                for ih:=1 to jb do vgrad[ih]:=xf[ih];


                f0:=fnev4(ev,th,prd);

                if f0>1e37 then
                begin
                     gf0:=f0;
                end;
                gf0:=f0;
                gprd:=ti0[4];
                gnak:=ti0[11];



                f2:=flam4(l,vgrad,ev,th,prd);
                l:=abs(l);
                svll:=l;
                while (f2>f0) and (l>1e-10) do
                begin
                     l:=l/10;
                     f2:=flam4(l,vgrad,ev,th,prd);
                end;
                bvht:=(l<1e-10);
                if (l<1e-10) {and false} {and( svll>1e-7)} or (abs(ev-1)<0.01) then
                begin


                     tsxf:=nsxf;
                     if bsxf then tsxf:=trs;
                     for ih:=1 to jb do
                     begin
                          sort3(sxf[ih]^,tsxf);
                          ih2:=trunc(tsxf/2);
                          if ih2<1 then ih2:=1;
                          if ih2>trs then ih2:=trs;


                          svgrad[ih]:=sxf[ih]^[ih2];
                     end;
                     {svgrad[2]:=0;}
                     min:=1e+38;
                     lp:=1e-5;
                     if svll>lp then lp:=svll;
                     {irn:=0;
                     repeat}
                           ir:=0;
                           repeat
                           begin
                                ir:=ir+1;
                                for ih:=1 to jb do ved[ih]:=2*random-1;
                                for ih:=1 to jb do
                                begin
                                     vgrad[ih]:=svgrad[ih]*ved[ih];
                                     {vgrad[ih]:=xf[ih]*ved[ih];}
                                end;
                                {fd8:=flam(lp,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia)-f0;}
                                fd8:=flam4(lp,vgrad,ev,th,prd)-f0;
                                if fd8<min then
                                begin
                                     min:=fd8;
                                     for ih:=1 to jb do mnved[ih]:=ved[ih];
                                end;
                           end;
                           until ((min<0) and (ir>chi)) or (ir>chi*10);
                           {if min>=0 then
                           begin
                                lp:=lp/10;
                           end;}
                           {irn:=irn+1;
                     until (min<0) or (irn>10);}
                     if min<0 then
                     begin
                          mnt:=1;
                          for ih:=1 to jb do
                          begin
                               xf[ih]:=svgrad[ih]*mnved[ih];
                               {xf[ih]:=xf[ih]*mnved[ih];}

                          end;

                          for ih:=1 to jb do vgrad[ih]:=xf[ih];
                          fd8:=min;

                          repeat
                                fd7:=fd8;
                                mnt:=mnt*2;

                                fd8:=flam4(lp*mnt,vgrad,ev,th,prd)-f0;
                          until (fd7<fd8) or (mnt*lp>1);
                          mnt:=mnt/2;
                          l:=lp*mnt;
                          vht:=vht+1;
                     end;

                     fd1:=0;
                end;
                if not(bvht) then
                begin
                     {f3:=flam(l/1.1,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);}
                     f3:=flam4(l/1.1,vgrad,ev,th,prd);
                     {f4:=flam(l*1.03,vgrad,iro0,iteta0,rot,tetat,r,vr,massa,pia);}
                     f4:=flam4(l*1.03,vgrad,ev,th,prd);
                     if (f3<f2) and (f3<f0) then l:=l/1.1;
                     if (f4<=f2) and (f4<f0) then l:=l*1.03;
                end;
                {l:=1e-2;}
                {if fd>1000 then l:=1e-1;}


                     jb:=0;
                     vyh:=true;

                     lp:=l;
                     if uev then
                     begin
                          jb:=jb+1;
                          ev:=ev-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<eps2);
                     end;
                     if uth then
                     begin

                          jb:=jb+1;
                          th:=th-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<eps2);
                     end;
                     {prd,ai,omm,omb}
                     if uprd then
                     begin
                          jb:=jb+1;
                          prd:=prd-lp*xf[jb];
                          vyh:=vyh and (abs(xf[jb])<eps2);

                     end;
                     if ev<0 then
                     begin
                          fd:=fd+1;
                     end;
                fd:=fd+1;
                {if massa<0 then
                begin
                     bvr:=sbvr;
                     bmassa:=sbmassa;
                     bpia:=sbpia;
                     pia:=spia;
                     gf0:=1e39;

                     exit;

                end;}

                GetTime(qHour, qMinute, qSecond, qSec100);
                qSecond:=qSec100 div 20;

                if (qSecond<>qqSecond) and not(bvnevyd) then
                begin

                {if (fd mod 1 =0) and not(bvnevyd) then
                begin}
                {abp,ev,th,prd,ai,omm,omb}

                gotoxy(1,1);
                if model then writeln('Model') else
                writeln;
                if uev  then pst:='     var'
                        else pst:=' not var'  ;
                writeln('e=',ev:10:5,pst);
                if uth  then pst:='     var'
                        else pst:=' not var'  ;
                writeln('T=',th+t0 :10:5,pst);
                if uprd then pst:='     var'
                        else pst:=' not var'  ;
                writeln('P=',prd  :10:5,pst);

                writeln('itr=',fd:5,' N.Eq.=',i:5,' d.l.=',vht:2 );

                writeln('RMS=',f0:15:20 );
                writeln('lamda=',lp:15 );
                {writeln('e=',ev:10:5 );}
                tkmjd:=tekmjdr-bmjd;
                tkmjd:=tkmjd*60*24;
                writeln('Minutes=',tkmjd:10:5 );
                writeln('N.tests=',iper);
                if iper=-10 then
                begin
                     iper:=-10;
                end;
                {if fixoi then
                begin
                     writeln('fixoi=true');
                end else
                begin
                     writeln('fixoi=false');
                end;}
                writeln(swds,' ',comp10);

                end;
                qqSecond:=qSecond;



                if not(uge) then
                begin
                     ku2:=ku2+1;
                     uge:=true;
                end;

                {if biskl then wvh;}
             ch:='p';

            if keypressed then ch:=readkey;
            if ch='a' then halt(0);
            tkmjd:=tekmjdr-bmjd;
            tkmjd:=tkmjd*60*24;



            nevyh:=(abs(ev-1.0)<0.001) or
                   (ev>1) or
                   {(r>3000) or}
                   (fd>t50000) or
                   {(abs(tetat)>100) or}
                   {((fd>200) and ((vht*10)>fd)) or
                   ((fd>10) and ((vht*2)>fd)) or}
                   {(massa>1000) or
                   (massa<0) or
                   ((abs(pia-1)<0.001) and (abs(vr)>10)) or
                   (abs(rot)>10) or}
                   (l<1e-10) or
                   (tkmjd>predmin) or
                   ((ti0[4]<1) and not((swds='14575-2125') and (comp10='Ba,Bb')))
                    ;
            if koi then
            begin
                   {nevyh:=nevyh or ((fd>200) and ((vht*10)>fd)) or
                   ((fd>10) and ((vht*2)>fd));}


            end;



          until (ch='a') or vyh or nevyh;

          if {(l>1e-10) and} (ev<50) then
          begin

          if biskl then
          begin
               maxiskl4(ev,th,prd,oevr,oevt,maxnev,mi);
               {(ev,th,prd:extended):extended;}
               ew:=false;
               begin
                    ew:=maxnev>porog;
                   if ((swds='06003-3102') and (comp10='CE   ')) or
                      ((swds='19072+2053')) or
                      ((swds='12283-6146') and (comp10='BC   ')) or
                      ((swds='06579-4417')) or
                      ((swds='13473+1727')) or
                      ((swds='03236-4005')) or
                      ((swds='08563-3707')) or
                      ((swds='16278-0822')) or
                      ((swds='14503+2355') and (comp10='BC   ')) or
                      ((swds='15041-0653')) or
                      ((swds='06032+5813')) or
                      ((swds='23075+2108')) or
                      ((swds='08122+1739') and (comp10='Ca,Cb'))
                   then
                   begin
                        ew:=maxnev>2.9;
                   end;


                    if ew then
                    begin
                         cvb2(mi,3);
                         {if nevyh then
                         begin
                               iro0:=siro0;
                               iteta0:=siteta0;
                               rot:=srot;
                               tetat:=stetat;
                               r:=sr;
                               massa:=svmassa;
                               vr:=svvr;

                         end;}
                         l:=ln;

                    end;
               end;
          end
          else ew:=false;
          end;
          wvh;


          until  (ch='a') or not(ew) or (tkmjd>predmin) ;

          {if biskl then vydres2(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt,oevr,oevt);}
          end;
          uci:=true;
          if ch='a' then uci:=false;
          {if l<1e-10 then uci:=false;}
          if ev>2 then uci:=false;
          if (r>3000) or (fd>t50000) or (abs(tetat)>1000) then uci:=false;
          if nevyh then uci:=false;

          uci:=true;

          if uci then
          begin
                sgmnk2(al,cv,xf,sgm,jb,i*2);

          end else
          begin
               gf0:=1e39;
                               iro0:=siro0;
                               iteta0:=siteta0;
                               rot:=srot;
                               tetat:=stetat;
                               r:=sr;
                               massa:=svmassa;
                               vr:=svvr;


          end;
      {massa:=svmassa;
      vr:=svvr;}

      for i:=1 to 8 do dispose(sxf[i]);
      {btvht:=((fd>200) and ((vht*10)>fd)) or ((fd>10) and ((vht*2)>fd));}




end;

procedure ishciriteracii(biskl:boolean);
var i,ih:integer;
    lp,rog2,ad,bd,h,rog0,tetag0:extended;
    rsx:text;
    ti0,tiro,titeta,tirot,titetat,tir,tivr,tim,tip,vgrad:das;
    tiro0,titeta0,tirot0,titetat0,tir0,tivr0,tim0,tip0:das;

    ch:char;
    f0,f1,f2,f3,f4,maxnev,ff1,ff2:extended;
    mi,itr:integer;
    sbvr,sbmassa,sbpia:boolean;
    spia,svmassa,svvr,raz,ln:extended;
    uge:boolean;
    siro0,siteta0,srot,stetat,sr,eps,eps2:extended;




begin
      bklpvd:=true;

      nevyh:=false;
      bzc:=true;

      btvht:=false;

      eps:=1e-7;
      eps2:=1e-6;
      l:=le;

      siro0:=iro0;
      siteta0:=iteta0;
      srot:=rot;
      stetat:=tetat;
      sr:=r;

      uci:=false;
      ln:=l;
      uge:=false;
      mbcir:=true;

      svmassa:=massa;
      svvr:=vr;
      sbvr:=bvr;
      sbmassa:=bmassa;
      sbpia:=bpia;
      spia:=pia;

      bvr:=false;
      bmassa:=false;
      bpia:=false;
      pia:=1;

      {br:=false;
      r:=sros*1000;}

      dro:=eps*180/pi;
      ch:='p';

         fd:=0;

          begin
                {massa:=massa-0.1;}
          repeat
          t3:=nt;


          repeat
                i:=0;
                mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia,ti0);
                ev:=ti0[2];
                {if ti0[2]>1 then
                begin
                      exit;
                end;}
                {if r<0.0001 then exit;}
                {if fd>550 then exit;}

        if bro0 then mpdq.cirtpvd(iro0-0.5*eps,iteta0,rot,tetat,r,vr,massa,et^[1],pia,tiro0);
        if bteta0 then mpdq.cirtpvd(iro0,iteta0-0.5*eps*180/pi,rot,tetat,r,vr,massa,et^[1],pia,titeta0);
        if brot then mpdq.cirtpvd(iro0,iteta0,rot-0.5*eps*0.01,tetat,r,vr,massa,et^[1],pia,tirot0);
        if btetat then mpdq.cirtpvd(iro0,iteta0,rot,tetat-0.5*eps,r,vr,massa,et^[1],pia,titetat0);
        if br then mpdq.cirtpvd(iro0,iteta0,rot,tetat,r-0.5*eps,vr,massa,et^[1],pia,tir0);
        {if bvr then mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,(vr-0.5*eps),massa,et[1],pia,tivr0);
        if bmassa then mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,vr,massa-0.5*eps,et[1],pia,tim0);}
        if bpia then  mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia-0.5*eps,tip0);

        if bro0 then mpdq.cirtpvd(iro0+0.5*eps,iteta0,rot,tetat,r,vr,massa,et^[1],pia,tiro);
        if bteta0 then mpdq.cirtpvd(iro0,iteta0+0.5*eps*180/pi,rot,tetat,r,vr,massa,et^[1],pia,titeta);
        if brot then mpdq.cirtpvd(iro0,iteta0,rot+0.5*eps*0.01,tetat,r,vr,massa,et^[1],pia,tirot);
        if btetat then mpdq.cirtpvd(iro0,iteta0,rot,tetat+0.5*eps,r,vr,massa,et^[1],pia,titetat);
        if br then mpdq.cirtpvd(iro0,iteta0,rot,tetat,r+0.5*eps,vr,massa,et^[1],pia,tir);
        {if bvr then mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,vr+0.5*eps,massa,et[1],pia,tivr);
        if bmassa then mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,vr,massa+0.5*eps,et[1],pia,tim);}
        if bpia then  mpdq.cirtpvd(iro0,iteta0,rot,tetat,r,vr,massa,et^[1],pia+0.5*eps,tip);
                i:=0;
                for id:=1 to t3 do
                begin
                     if gbop^[id]=0 then
                     begin
                          tcir2(ti0,tiro0,tiro,titeta0,titeta,
                               tirot,tirot0,titetat0,titetat,
                                tir0,tir,tivr0,tivr,tim0,tim,tip0,tip,
                                rog0,tetag0,eps,raz,i,
                                jb,id,
                                et,ro,teta,
                                rog,tetag,rog1,tetag1,dro,
                                bro0,bteta0,brot,btetat,br,bvr,bmassa,bpia,
                                al
                                 );
                     end;
                end;
                {halt(0);}
                {close(rsx);}
                ptlst4(al,xf,cv,i*2,jb);
                {sgmnk(a,cv,xf,sgm,jb,i*2);}

                for ih:=1 to jb do vgrad[ih]:=xf[ih];

                f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
                gf0:=f0;
                gprd:=ti0[4];
                gnak:=ti0[11];
                {if fd=348 then
                begin
                     ih:=0;
                end;}


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
                if massa<0 then
                begin
                     bvr:=sbvr;
                     bmassa:=sbmassa;
                     bpia:=sbpia;
                     pia:=spia;
                     mbcir:=false;
                     gf0:=1e39;

                     exit;

                end;



                {if fd mod 1 =0 then}
                if (fd mod 1 =0) and not(bvnevyd) then

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
                {readln;}
                if not(uge) then
                begin
                     ku2:=ku2+1;
                     uge:=true;
                end;

                if biskl then wvh;
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

            if (l<1e-10) or (ev>2) or (abs(tetat)>20) then
            begin
                 {vyhtup(iro0,iteta0,rot,tetat,r,vr,massa,pia);}
                  vyh:=true;
                 l:=1e-15;
            end;
            if ch='a' then halt(0);
            nevyh:=f0>0.9e38;


          until (ch='a') or vyh or (r>3000) or (fd>t50000) or (abs(tetat)>1000) or nevyh;
          {until vyh;}
          {vydres(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt);;;;}

          {nevotap(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt);;;;}
          if {(l>1e-10) and} (ev<50) then
          begin

          if biskl then
          begin
               {maxnev2(iro0,iteta0,rot,tetat,r,vr,massa,pia,oevr,oevt);
               iskl2(iro0,iteta0,rot,tetat,r,vr,massa,pia,oevr,oevt,maxnev,mi);}
               maxiskl(iro0,iteta0,rot,tetat,r,vr,massa,pia,oevr,oevt,maxnev,mi);
               ew:=false;
               if (mi>1) and (mi<nt) then
               begin
                    ew:=maxnev>porog;
                    if ew then
                    begin
                         cvb2(mi,2);
                         if l<1e-10 then
                         begin
                               iro0:=siro0;
                               iteta0:=siteta0;
                               rot:=srot;
                               tetat:=stetat;
                               r:=sr;
                         end;
                         {l:=ln;}

                    end;
               end;
          end
          else ew:=false;
          {ew:=false;}
          {if biskl then wvh;}
          end;
          wvh;
          {if ((l<1e-10) and ew) or (ev>50) then ch:='a';}


          until  (ch='a') or not(ew) or nevyh;
          {if biskl then vydres2(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt,oevr,oevt);}
          end;
          uci:=true;
          if nevyh then uci:=false;
          if ch='a' then uci:=false;
          {if l<1e-10 then uci:=false;}
          if ev>50 then uci:=false;
          if (r>2500) or (fd>t50000) or (abs(tetat)>1000) then uci:=false;
          tekr:=r;
          {if nevyh then uci:=false;}
          if uci then
          begin
                sgmnk2(al,cv,xf,sgm,jb,i*2);

          end else
          begin
               gf0:=1e39;
                               iro0:=siro0;
                               iteta0:=siteta0;
                               rot:=srot;
                               tetat:=stetat;
                               r:=sr;

          end;
      {massa:=svmassa;
      vr:=svvr;}

      bvr:=sbvr;
      bmassa:=sbmassa;
      bpia:=sbpia;
      pia:=spia;
      mbcir:=false;



end;
procedure ciriteracii(biskl:boolean);
var i,ih:integer;
begin
     bklpvd:=true;

     br:=true;
     ishciriteracii(biskl);
     if tekr>2500 then
     begin
          r:=1000;
          br:=false;
          ishciriteracii(biskl);
          br:=true;
     end;
end;



begin
end.
