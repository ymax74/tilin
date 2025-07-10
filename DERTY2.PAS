unit derty2;

interface
         uses sbop,mpdq,ubrk,mnk,mpvd,crt,derty;
         var
         xfk:array[1..12] of extended;
         {procedure odflM3d(pst:string);}
         procedure odflM4d(pst:string);
         procedure odflM5d(pst:string);
         {procedure vychosh(pst:string);}
         procedure itgtr;


implementation
const botl=false;
procedure odflM3d(pst:string);
const t3=3;
var id,jd,kd,kt,i,mi,mi2,mi3,mi4,nt2,j,t10,mi5:integer;
min,max:array[1..t3] of extended;
bmx:boolean;
mix,mix2:array[1..t3] of integer;
rop,st0,mnt,mxt,mxr,raz,raz2,xc,yc,mnt2,svvr,th2:extended;
ti0:das;
t2,diro0,diteta0,drot,dtetat,dr,dvr,dmassa,oev,mint,mprd,mev,mth,mabp,momm,momb,mai:extended;
rs,rs2:text;
aa2,bb2,ff2,gg2,comb,somb,comm,somm,ci,f0,tt100,f1,tkmjd:extended;
bcm:boolean;
cai,tai,clomb,tomb:extended;
asn1,asn2,asn3,asn4,asn5:ats;

begin


      isuge:=false;
      tfl:=pst;
      niskl:=0;
      ClrScr;

      mi:=-1;
      new(asn1);new(asn2);new(asn3);new(asn4);new(asn5);
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
      if (nt>10) and (mi>0) and bcm then
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
            mi5:=-1;
             tkmjd:=tekmjdr-bmjd;
             tkmjd:=tkmjd*60*24;
             writeln('tkmjd=',tkmjd:10:5 );


            {if (pst12='S01350S2.955') or (pst12='S17066N0.039') then tt100:=10;}
            for i:=1 to t16383 do
            begin
                 ev:=random;
                 {ev:=0.784;}
                 prd:=random*tt100;
                 th:=random*prd-prd*0.5;
                 GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
                 GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
                 {ai:=atai[mi];
                 omb:=atomb[mi];}
                 asn1^[i]:=oev;
                 asn2^[i]:=ev;
                 asn3^[i]:=prd;
                 asn4^[i]:=th;
                 asn5^[i]:=GetUgol(abp,ev,th,prd,ai,omb,omm,
                                   atabp[mi],atev[mi],atth[mi],atprd[mi],atai[mi],atomm[mi],atomb[mi]);
                 if asn5^[i]>pi/2 then asn5^[i]:=pi-asn5^[i];

                 if (oev<mint)  then
                 begin

                      mint:=oev;
                      mprd:=prd;
                      mev:=ev;
                      mth:=th;
                      mi5:=i;
                 end;
           end;
           sort7(asn1,asn2,asn3,asn4,asn5,i);
           assign(rs,'tsmnt2.res');
           rewrite(rs);
           t10:=10;
           for i:=1 to t10 do
           begin
                ev:=asn2^[i];
                prd:=asn3^[i];
                th:=asn4^[i];
                {elemiter2(false);}
                GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
                GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
                GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);
                vr:=GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia);

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
           elemiter2(true);



           GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
           GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
           GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);
           vr:=GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia);
           f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);

           mbcir:=false;
           vydnev(1);
           i:=0;
           for id:=1 to nt do
           begin
                if gbop^[id]<>0 then i:=i+1;
           end;

           writeln(rstr,swds,' ',comp10,' ',abp:20:15,' ',ev:20:15,' ',th+t0:20:15,' ',
                         prd:20:15,' ',ai:20:15,' ',omm:20:15,' ',omb:20:15,' ',
                         t0:20:15,' ',
                         iro0:20:15,' ',iteta0:20:15,' ',rot:20:15,' ',
                         tetat:20:15,' ',r:20:15,' ',vr:20:15,' ',massa:20:15,' ',f0:20:15,' ',nt:4,' ',i:4);
           flush(rstr);




           close(rs);

           vydres4;


           sort7(asn5,asn1,asn2,asn3,asn4,i);

           sort7(asn1,asn2,asn3,asn4,asn5,trunc(i/16));
             tkmjd:=tekmjdr-bmjd;
             tkmjd:=tkmjd*60*24;
             writeln('tkmjd=',tkmjd:10:5 );


           {sort7(asn1,asn2,asn3,asn4,asn5,i);}
           {writeln('ku');
           readln;}
           assign(rs,'tsmnt2.res');
           rewrite(rs);
           mint:=1e39;
           t10:=10;
           {if pst12='S08122N1.739' then t10:=1;
           t10:=1;}

           for i:=1 to t10 do
           begin
                ev:=asn2^[i];
                prd:=asn3^[i];
                th:=asn4^[i];
                {elemiter2(false);}
                GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
                GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
                GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);
                vr:=GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia);

                uomb:=false;
                uai:=false;

                elemiter(false);

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




                {readln;}
                writeln(rs,i:3,' ',gf0:8:5,' ',ev:8:5,' ',prd:15:5,' ',th+t0:15:5);
                 if (gf0<mint)  then
                 begin
                      mint:=gf0;
                      mprd:=prd;
                      mev:=ev;
                      mth:=th;
                      mabp:=abp;
                      momm:=omm;
                      momb:=omb;
                      mai:=ai;

                 end;

           end;
           prd:=mprd;
           ev:=mev;
           th:=mth;
           abp:=mabp;
           omm:=momm;
           omb:=momb;
           ai:=mai;
           GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);
           vr:=GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia);
           f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
               i:=0;
               for id:=1 to nt do
               begin
                    if gbop^[id]<>0 then i:=i+1;
               end;



           writeln(rstr2,swds,' ',comp10,' ',abp:20:15,' ',ev:20:15,' ',th+t0:20:15,' ',
                         prd:20:15,' ',ai:20:15,' ',omm:20:15,' ',omb:20:15,' ',
                         t0:20:15,' ',
                         iro0:20:15,' ',iteta0:20:15,' ',rot:20:15,' ',
                         tetat:20:15,' ',r:20:15,' ',vr:20:15,' ',massa:20:15,' ',f0:20:15,' ',nt:4,' ',i:4);
           flush(rstr2);

           close(rs);
           vydres4;


      end;
      dispose(asn1);
      dispose(asn2);
      dispose(asn3);
      dispose(asn4);
      dispose(asn5);



end;
procedure mappend(var rs:text);
begin
{$I-}
     append(rs);
     if ioresult<>0 then rewrite(rs);
{$I+}
end;
function GetRanPer2:extended;
var rprd,mn,nev,prd:extended;
t10,j:integer;
begin
          rprd:=random*500;
      { rprd - случайное число от 0 до 500}
          t10:=11;
          nev:=0;
          mn:=1;
          for j:=0 to t10 do
          begin
               nev:=nev+mn*xfk[j+1];
               mn:=mn*rprd;
          end;
          GetRanPer2:=exp(nev);

end;

function GetRanPer:extended;
var rprd,mn,nev,prd:extended;
t10,j:integer;
begin
          rprd:=GetRanPer2;
          if (swds='04400+5328') and (comp10='AB   ') then
          begin
                while rprd>50 do
                begin
                     rprd:=GetRanPer2;
                end;
          end;
          if (swds='22038+6438') and (comp10='Aa,Ab') then
          begin
                while rprd>5 do
                begin
                     rprd:=GetRanPer2;
                end;
          end;
          if (swds='14575-2125') and (comp10='Ba,Bb') then
          begin
                while rprd>3 do
                begin
                     rprd:=GetRanPer2;
                end;
          end;
          GetRanPer:=rprd;

end;

procedure odflM4d(pst:string);
const t3=3;
var id,jd,kd,kt,i,mi,mi2,mi3,mi4,nt2,j,t10,mi5,mi6:integer;
min,max:array[1..t3] of extended;
bmx:boolean;
mix,mix2:array[1..t3] of integer;
rop,st0,mnt,mxt,mxr,raz,raz2,xc,yc,mnt2,svvr,th2:extended;
ti0:das;
t2,diro0,diteta0,drot,dtetat,dr,dvr,dmassa,oev,mint,mprd,mev,mth,mabp,momm,momb,mai:extended;
rs,rs2:text;
aa2,bb2,ff2,gg2,comb,somb,comm,somm,ci,f0,tt100,f1,tkmjd,pi2:extended;
bcm:boolean;
cai,tai,clomb,tomb:extended;
asn1,asn2,asn3,asn4,asn5:ats;

begin


      isuge:=false;
      tfl:=pst;
      niskl:=0;
      ClrScr;

      mi:=-1;
      new(asn1);new(asn2);new(asn3);new(asn4);new(asn5);
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
      mi6:=-1;
      FOR I:=1 TO 15 do
      begin
           if swds=atswds[i] then mi6:=i;
      end;
      {if btesorb then}
      begin
           mi:=1;
           mi6:=1;
           bcm:=true;
      end;


      if (nt>3) and (mi>0) and bcm and (mi6>0) then
      begin


           {clan(m1,m2,teta,x,y,ro,nnt,ain,aex,aex2,gbop,gbop2,gbop3,aex3);}
           if not(fixoi2) and not(bVneIskl) then
           begin
                for id:=1 to nt do
                begin
                     gbop^[id]:=0;
                     gbop3^[id]:=0;
                end;
           end{ else
           begin
                 cai:=atai[mi6];
                 omb:=atomb[mi6];
                 pi2:=pi/2;
                 if ((cai>pi2) and (ai<pi2)) or ((cai<pi2) and (ai>pi2)) then
                 begin
                      cai:=pi-cai;
                 end;
                 ai:=cai;
           end};
           if (swds='06173+0506') and (comp10='Aa,Ab') then
           begin
                gbop^[5]:=1;
           end;
           if ((swds='19072+2053'))
           then
           begin
                for id:=1 to nt do
                begin
                     IF Aref^[id]='USN1974 ' then
                     begin
                          gbop^[id]:=1;
                     end;
                end;

           end;
           if ((swds='14575-2125'))
           then
           begin
                for id:=1 to nt do
                begin
                     IF Aref^[id]='Prr1990 ' then
                     begin
                          gbop^[id]:=1;
                     end;
                end;

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

            uth:=true;
            mint:=1e39;
            ev:=j*0.05;
            if ev>=1 then ev:=0.999;
            tt100:=2000;
            if btesorb then tt100:=70;
            mi5:=-1;
             tkmjd:=tekmjdr-bmjd;
             tkmjd:=tkmjd*60*24;
             writeln('tkmjd=',tkmjd:10:5 );


            {if (pst12='S01350S2.955') or (pst12='S17066N0.039') then tt100:=10;}
            for i:=1 to t16383 do
            begin
                 ev:=random;
                 if ev<0.99 then
                 begin
                 prd:=GetRanPer;
                 {prd:=random*tt100;}
                 th:=random*prd-prd*0.5;
                 GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
                 GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
                 {ai:=atai[mi];
                 omb:=atomb[mi];}
                 asn1^[i]:=oev;
                 asn2^[i]:=ev;
                 asn3^[i]:=prd;
                 asn4^[i]:=th;
                 asn5^[i]:=0;

                 {asn5^[i]:=GetUgol(abp,ev,th,prd,ai,omb,omm,
                                   atabp[mi],atev[mi],atth[mi],atprd[mi],atai[mi],atomm[mi],atomb[mi]);
                 aa2:=GetUgol(abp,ev,th,prd,ai,omb+pi,omm+pi,
                                   atabp[mi],atev[mi],atth[mi],atprd[mi],atai[mi],atomm[mi],atomb[mi]);
                 if asn5^[i]>pi/2 then
                 begin
                      asn5^[i]:=pi-asn5^[i];
                 end;}

                 if (oev<mint)  then
                 begin

                      mint:=oev;
                      mprd:=prd;
                      mev:=ev;
                      mth:=th;
                      mi5:=i;
                 end;
                 end else
                 begin
                        asn1^[i]:=1e5;
                        asn2^[i]:=ev;
                        asn3^[i]:=prd;
                        asn4^[i]:=th;
                        asn5^[i]:=0;


                 end;
           end;
           sort7(asn1,asn2,asn3,asn4,asn5,i);
           p1:=pos('.',mflnm);
           pst:=copy(mflnm,1,p1)+'nt2';

           {assign(rs,'tsmnt2.res');}
           {assign(rs,pst);
           mappend(rs);}
           t10:=10;
           if botl then t10:=1;
           i:=0;
           {for i:=1 to t10 do}
           {writeln(rs);}
           mint:=1e39;
           {i:=iper;}
           if bitrves then t10:=5;
           while (i<t10) and (tkmjd<predmin) do
           begin
                i:=i+1;
                ev:=asn2^[i];
                prd:=asn3^[i];
                th:=asn4^[i];
                {elemiter2(false);}
                GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
                GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
                GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);
                vr:=GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia);
                iper:=i;
                if fixoi2 then iper:=t10+i;


                elemiter2(false);
                {elemiter2(true);}
                {readln;}
                 tkmjd:=tekmjdr-bmjd;
                 tkmjd:=tkmjd*60*24;

                iper:=0;
                if fixoi then iper:=1;
                {writeln(rs,swds,' ',comp10,' ',iper,' ',i:3,' ',gf0:8:5,' ',
                ev:8:5,' ',prd:15:5,' ',th+t0:15:5,' ',tkmjd:8:2);}
                 if (gf0<mint)  then
                 begin
                      mint:=gf0;
                      mprd:=prd;
                      mev:=ev;
                      mth:=th;
                 end;


           end;
           {close(rs);}
           bmjd:=tekmjdr;

           prd:=mprd;
           ev:=mev;
           th:=mth;
           iper:=-i;
           if not(fixoi2) then
           begin
                 {elemiter2(true);}
                 {elemiter2(false);}
                 elemiter2(not(bVneIskl));
           end else
           begin
                 elemiter2(false);
           end;



           GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
           GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
           glev:=ev;
           glth:=th;
           glprd:=prd;
           glabp:=abp;
           glai:=ai;
           glomb:=omb;
           glomm:=omm;
           GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);
           vr:=GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia);
           f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);

           mbcir:=false;
           {if not(fixoi) then
           begin
                vydnev(1);
           end;}
           i:=0;
           for id:=1 to nt do
           begin
                if gbop^[id]<>0 then i:=i+1;
           end;
           {if not(fixoi) then
           begin

                 writeln(rstr,swds,' ',comp10,' ',abp:20:15,' ',ev:20:15,' ',th+t0:20:15,' ',
                         prd:20:15,' ',ai:20:15,' ',omm:20:15,' ',omb:20:15,' ',
                         t0:20:15,' ',
                         iro0:20:15,' ',iteta0:20:15,' ',rot:20:15,' ',
                         tetat:20:15,' ',r:20:15,' ',vr:20:15,' ',massa:20:15,' ',f0:20:15,' ',nt:4,' ',i:4);
                flush(rstr);
           end else
           begin
                 writeln(rstr2,swds,' ',comp10,' ',abp:20:15,' ',ev:20:15,' ',th+t0:20:15,' ',
                         prd:20:15,' ',ai:20:15,' ',omm:20:15,' ',omb:20:15,' ',
                         t0:20:15,' ',
                         iro0:20:15,' ',iteta0:20:15,' ',rot:20:15,' ',
                         tetat:20:15,' ',r:20:15,' ',vr:20:15,' ',massa:20:15,' ',f0:20:15,' ',nt:4,' ',i:4);
                flush(rstr2);


           end;}



      end;
      dispose(asn1);
      dispose(asn2);
      dispose(asn3);
      dispose(asn4);
      dispose(asn5);



end;

procedure odflM5d(pst:string);
const t3=3;
var id,jd,kd,kt,i,mi,mi2,mi3,mi4,nt2,j,t10,mi5,mi6:integer;
min,max:array[1..t3] of extended;
bmx:boolean;
mix,mix2:array[1..t3] of integer;
rop,st0,mnt,mxt,mxr,raz,raz2,xc,yc,mnt2,svvr,th2:extended;
ti0:das;
t2,diro0,diteta0,drot,dtetat,dr,dvr,dmassa,oev,mint,mprd,mev,mth,mabp,momm,momb,mai:extended;
rs,rs2:text;
aa2,bb2,ff2,gg2,comb,somb,comm,somm,ci,f0,tt100,f1,tkmjd,pi2,oev1,oev2,oev3,oev4,oev5,vr1,vr2,vr3:extended;
bcm:boolean;
cai,tai,clomb,tomb:extended;
asn1,asn2,asn3,asn4,asn5:ats;
miro0,miteta0,mrot,mtetat,mr,mvr,mmassa,mpia:extended;

begin


      isuge:=false;
      tfl:=pst;
      niskl:=0;
      ClrScr;

      mi:=-1;
      new(asn1);new(asn2);new(asn3);new(asn4);new(asn5);
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
      mi6:=-1;
      FOR I:=1 TO 15 do
      begin
           if swds=atswds[i] then mi6:=i;
      end;
      {if btesorb then}
      begin
           mi:=1;
           mi6:=1;
           bcm:=true;
      end;


      if (nt>3) and (mi>0) and bcm and (mi6>0) then
      begin


           {clan(m1,m2,teta,x,y,ro,nnt,ain,aex,aex2,gbop,gbop2,gbop3,aex3);}
           if not(fixoi2) then
           begin
                for id:=1 to nt do
                begin
                     gbop^[id]:=0;
                     gbop3^[id]:=0;
                end;
           end{ else
           begin
                 cai:=atai[mi6];
                 omb:=atomb[mi6];
                 pi2:=pi/2;
                 if ((cai>pi2) and (ai<pi2)) or ((cai<pi2) and (ai>pi2)) then
                 begin
                      cai:=pi-cai;
                 end;
                 ai:=cai;
           end};
           if (swds='06173+0506') and (comp10='Aa,Ab') then
           begin
                gbop^[5]:=1;
           end;
           if ((swds='19072+2053'))
           then
           begin
                for id:=1 to nt do
                begin
                     IF Aref^[id]='USN1974 ' then
                     begin
                          gbop^[id]:=1;
                     end;
                end;

           end;
           if ((swds='14575-2125'))
           then
           begin
                for id:=1 to nt do
                begin
                     IF Aref^[id]='Prr1990 ' then
                     begin
                          gbop^[id]:=1;
                     end;
                end;

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

            uth:=true;
            mint:=1e39;
            ev:=j*0.05;
            if ev>=1 then ev:=0.999;
            tt100:=2000;
            if btesorb then tt100:=70;
            mi5:=-1;
             tkmjd:=tekmjdr-bmjd;
             tkmjd:=tkmjd*60*24;
             writeln('tkmjd=',tkmjd:10:5 );


            {if (pst12='S01350S2.955') or (pst12='S17066N0.039') then tt100:=10;}
            for i:=1 to t16383 do
            begin
                 ev:=random;
                 if ev<0.99 then
                 begin
                 prd:=GetRanPer;
                 {prd:=random*tt100;}
                 th:=random*prd-prd*0.5;
                 GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
                 GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
                 GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);
                 {if i=5 then
                 begin
                      i:=5;
                 end;}

                 {oev1:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);
                 oev2:=fnev(iro0,iteta0,rot,tetat,r,-vr,massa,pia);}
                 {oev3:=fnev4(ev,th,prd);}
                 oev4:=fnev(iro0,iteta0,rot,tetat,r,ovr,omas,pia);
                 oev5:=fnev(iro0,iteta0,rot,tetat,r,-ovr,omas,pia);
                 if oev5<oev4 then oev4:=oev5;
                 {writeln(i);}


                 {ai:=atai[mi];
                 omb:=atomb[mi];}
                 asn1^[i]:=oev4;
                 {asn1^[i]:=oev3;}
                 asn2^[i]:=ev;
                 asn3^[i]:=prd;
                 asn4^[i]:=th;
                 asn5^[i]:=0;

                 {asn5^[i]:=GetUgol(abp,ev,th,prd,ai,omb,omm,
                                   atabp[mi],atev[mi],atth[mi],atprd[mi],atai[mi],atomm[mi],atomb[mi]);
                 aa2:=GetUgol(abp,ev,th,prd,ai,omb+pi,omm+pi,
                                   atabp[mi],atev[mi],atth[mi],atprd[mi],atai[mi],atomm[mi],atomb[mi]);
                 if asn5^[i]>pi/2 then
                 begin
                      asn5^[i]:=pi-asn5^[i];
                 end;}

                 if (oev<mint)  then
                 begin

                      mint:=oev4;
                      mprd:=prd;
                      mev:=ev;
                      mth:=th;
                      mi5:=i;
                 end;
                 end else
                 begin
                        asn1^[i]:=1e5;
                        asn2^[i]:=ev;
                        asn3^[i]:=prd;
                        asn4^[i]:=th;
                        asn5^[i]:=0;


                 end;
           end;
           sort7(asn1,asn2,asn3,asn4,asn5,i);
           p1:=pos('.',mflnm);
           pst:=copy(mflnm,1,p1)+'nt2';

           {assign(rs,'tsmnt2.res');}
           assign(rs,pst);
           mappend(rs);
           t10:=10;
           {t10:=1;}
           i:=0;
           {for i:=1 to t10 do}
           writeln(rs);
           mint:=1e39;
           {i:=iper;}
           if bitrves then t10:=5;
           {t10:=1;}
           while (i<t10) and (tkmjd<predmin) do
           begin
                i:=i+1;
                ev:=asn2^[i];
                prd:=asn3^[i];
                th:=asn4^[i];
                {elemiter2(false);}
                GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
                GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
                GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);
                massa:=omas;
                vr:=ovr;
                vr:=GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia);
                iper:=i;
                if fixoi2 then iper:=t10+i;
                bro0:=true;
                bteta0:=true;
                btetat:=true;
                brot:=true;
                br:=true;
                bvr:=false;
                bmassa:=false;

                bvr:=true;
                bmassa:=true;


                bpia:=false;
                riteracii(false);

                {elemiter2(false);}
                {elemiter2(true);}
                {readln;}
                 tkmjd:=tekmjdr-bmjd;
                 tkmjd:=tkmjd*60*24;

                iper:=0;
                if fixoi then iper:=1;
                writeln(rs,swds,' ',comp10,' ',iper,' ',i:3,' ',gf0:8:5,' ',ev:8:5,' ',prd:15:5,' ',th+t0:15:5,' ',tkmjd:8:2);
                 if (gf0<mint)  then
                 begin
                      mint:=gf0;
                      miro0:=iro0;
                      miteta0:=iteta0;
                      mrot:=rot;
                      mtetat:=tetat;
                      mr:=r;
                      mvr:=vr;
                      mmassa:=massa;
                      mpia:=pia;
                      mprd:=prd;
                      mev:=ev;
                      mth:=th;
                 end;


           end;
           close(rs);
           bmjd:=tekmjdr;

           prd:=mprd;
           ev:=mev;
           th:=mth;
                      iro0:=miro0;
                      iteta0:=miteta0;
                      rot:=mrot;
                      tetat:=mtetat;
                      r:=mr;
                      vr:=mvr;
                      massa:=mmassa;
                      pia:=mpia;

           iper:=-i;
           {if not(fixoi2) then
           begin
                 elemiter2(true);
           end else
           begin
                 elemiter2(false);
           end;}

           riteracii(false);
           bklpvd:=true;



           {GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
           GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
           glev:=ev;
           glth:=th;
           glprd:=prd;
           glabp:=abp;
           glai:=ai;
           glomb:=omb;
           glomm:=omm;
           GetPvdFrmEl(abp,ev,th,prd,ai,omm,omb,0,iro0,iteta0,rot,tetat,r,vr,massa);
           vr:=GetZnVr(iro0,iteta0,rot,tetat,r,vr,massa,pia);}
           f0:=fnev(iro0,iteta0,rot,tetat,r,vr,massa,pia);

           mbcir:=false;
           {if not(fixoi) then
           begin
                vydnev(1);
           end;}
           i:=0;
           for id:=1 to nt do
           begin
                if gbop^[id]<>0 then i:=i+1;
           end;
           {if not(fixoi) then
           begin

                 writeln(rstr,swds,' ',comp10,' ',abp:20:15,' ',ev:20:15,' ',th+t0:20:15,' ',
                         prd:20:15,' ',ai:20:15,' ',omm:20:15,' ',omb:20:15,' ',
                         t0:20:15,' ',
                         iro0:20:15,' ',iteta0:20:15,' ',rot:20:15,' ',
                         tetat:20:15,' ',r:20:15,' ',vr:20:15,' ',massa:20:15,' ',f0:20:15,' ',nt:4,' ',i:4);
                flush(rstr);
           end else
           begin
                 writeln(rstr2,swds,' ',comp10,' ',abp:20:15,' ',ev:20:15,' ',th+t0:20:15,' ',
                         prd:20:15,' ',ai:20:15,' ',omm:20:15,' ',omb:20:15,' ',
                         t0:20:15,' ',
                         iro0:20:15,' ',iteta0:20:15,' ',rot:20:15,' ',
                         tetat:20:15,' ',r:20:15,' ',vr:20:15,' ',massa:20:15,' ',f0:20:15,' ',nt:4,' ',i:4);
                flush(rstr2);


           end;}



      end;
      dispose(asn1);
      dispose(asn2);
      dispose(asn3);
      dispose(asn4);
      dispose(asn5);



end;
(*procedure vychosh(pst:string);
var i,j,k,mi5:integer;
ofl:string;
dt:text;
swds2,comp10a:string[10];
abpq,evq,thq,prdq,aiq,ommq,ombq,
t0q,iro0q,iteta0q,rotq,tetatq,rq,vrq,massaq,f0q,dro,dte,oev,tkmjd,tbk:extended;
ntq,iq,id:integer;
sref:string[48];
ti0:das;
rs:text;
asn1,asn2,asn3,asn4,asn5:ats;
mint,mprd,mev,mth:extended;
begin
     i:=pos('.',pst);
     tkmjd:=tekmjdr-bmjd;
     tkmjd:=tkmjd*60*24;
     new(asn1);new(asn2);new(asn3);new(asn4);new(asn5);


     begin
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
            abpq:=glabp;evq:=glev;thq:=glth;prdq:=glprd;
            aiq:=glai;ommq:=glomm;ombq:=glomb; f0q:=glf0q;

               {getelem(abpq,evq,thq,prdq,aiq,ommq,ombq,et^[1],pia,ti0);}
               for id:=1 to nt do
               begin
                   {if gbop^[id]=0 then}
                   begin
                         mpdq.rtti(ti0,et^[id],rog,tetag);
                         dro:=rog-roo^[id];
                         dte:=(tetag-tetao^[id]);
                         ogr2(dte);
                         roo^[id]:=rog;
                         tetao^[id]:=tetag;
                         k:=0;
                   end;
               end;
               assign(rs,'model.orb');
               rewrite(rs);
               write(rs,
'           a                   e                T                              P                    i                omega');
              writeln(rs,'               Node                  RMS   ');
{   0.509402078780292    0.156459913618909 2000.706118050688280          15.597608205761859   43.815508557203957  107.73769
 5728460505  173.528158033733258    0.012245400204052}
               k:=0;
               while (k<200) and (tkmjd<24*60) do
               begin
                    k:=k+1;

                    for id:=1 to nt do
                    begin
                         {if gbop^[id]=0 then}
                         begin
                              rog:=ranorm(f0q);
                              ro^[id]:=roo^[id]+rog;
                              rog:=ranorm(f0q)/roo^[id];
                              rog:=rog*180/pi;
                              teta^[id]:=tetao^[id]+rog;
                              x^[id]:=ro^[id]*sin(teta^[id]*pi/180);
                              y^[id]:=ro^[id]*cos(teta^[id]*pi/180);


                         end;
                    end;
                    ev:=evq;th:=thq;prd:=prdq;
                    abp:=abpq;;ai:=aiq;omm:=ommq;omb:=ombq;
                    uabp:=true;
                    uev:=true;
                    uth:=true;
                    uprd:=true;
                    uai:=true;
                    uomm:=true;
                    uomb:=true;
                    iper:=k;
                    mint:=1e39;
                    for i:=1 to t16383 do
                    begin
                         ev:=random;
                         if ev<0.99 then
                         begin
                              prd:=GetRanPer;
                              {prd:=random*tt100;}
                              th:=random*prd-prd*0.5;
                              GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
                              GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
                              {ai:=atai[mi];
                              omb:=atomb[mi];}
                              asn1^[i]:=oev;
                              asn2^[i]:=ev;
                              asn3^[i]:=prd;
                              asn4^[i]:=th;
                              asn5^[i]:=0;
                              if (oev<mint)  then
                              begin

                                    mint:=oev;
                                    mprd:=prd;
                                    mev:=ev;
                                    mth:=th;
                                    mi5:=i;
                              end else
                              begin
                                   asn1^[i]:=1e5;
                                   asn2^[i]:=ev;
                                   asn3^[i]:=prd;
                                   asn4^[i]:=th;
                                   asn5^[i]:=0;


                              end;
                         end;
                    end;
                    model:=true;
                    ev:=mev;
                    th:=mth;
                    prd:=mprd;

                    elemiter2(false);
                    model:=false;

                    id:=0;
                    GetPrmEl(ev,th,prd,aa,bb,ff,gg,oev);
                    GetElFrTi(abp,ai,omb,omm,aa,bb,ff,gg);
                    tbk:=180-ai*180/pi;
                    if omb>pi then
                    begin
                          omb:=omb-pi;
                          omm:=omm-pi;
                          if omm<0 then omm:=2*pi+omm;
                    end;
                    tbk:=180-ai*180/pi;




                    writeln(rs,abp:20:15,' ',ev:20:15,' ',th+t0:20:15,' ',
                    prd:27:15,' ',tbk:20:15,' ',omm*180/pi:20:15,' ',omb*180/pi:20:15,' ',
                         gf0:20:15);
                    flush(rs);
                    tkmjd:=tekmjdr-bmjd;
                    tkmjd:=tkmjd*60*24;


               end;
               close(rs);

          i:=0;
     end;
     dispose(asn1);
     dispose(asn2);
     dispose(asn3);
     dispose(asn4);
     dispose(asn5);

end;*)
procedure itgtr;
var srd,oev90:extended;
ofl:string;
abpq,evq,thq,prdq,aiq,ommq,ombq,
t0q,iro0q,iteta0q,rotq,tetatq,rq,vrq,massaq,f0q,dro,dte,oev,tkmjd:extended;
ntq,iq,id:integer;

abpw,evw,thw,prdw,aiw,ommw,ombw,
t0w,iro0w,iteta0w,rotw,tetatw,rw,vrw,massaw,f0w:extended;
ntw,iw:integer;

abpy,evy,thy,prdy,aiy,ommy,omby,
t0y,iro0y,iteta0y,roty,tetaty,ry,vry,massay,f0y:extended;
nty,iy:integer;

swds2,comp10a,compy:string[10];
rs,rs2:text;
ug1,ug2:extended;
wcomp:string[5];

begin
      assign(rs,'trotv.res');
      rewrite(rs);
      assign(rs2,'trorb.res');
      rewrite(rs2);
      spfl(ovds+'\as91\a??????k.res',sf,ka);

      for ia:=1 to ka do
      begin
           assign(dt,ovds+'\as91\'+sf^[ia]);
           reset(dt);
           km:=0;
           while not(eof(dt)) do
           begin
                km:=km+1;
                readln(dt,abp,ev,th,prd,ai,omm,omb,
                      aex^[km]);

           end;
           sort(aex,km);
           writeln(sf^[ia],' ',km);
           if sf^[ia]='A100035K.RES' then
           begin
                j:=0;
           end;
           if km mod 2 = 0 then
           begin
                 i:=trunc(km/2);
                 srd:=(aex^[i]+aex^[i+1])/2;
           end else
           begin
                 i:=trunc(km/2)+1;
                 srd:=aex^[i];
           end;
           i:=trunc(km*0.95+0.5);
           oev90:=aex^[i];
           oev90:=oev90-srd;
           close(dt);

           i:=pos('.',sf^[ia]);
           ofl:=copy(sf^[ia],1,i-2);
           assign(dt,ovds+'\as91\'+ofl+'.txt');
           {writeln(ofl);}
           reset(dt);
           read(dt,pst12);
           read(dt,ch,comp10);
           swds:=copy(pst12,2,7)+copy(pst12,10,3);
           if swds[6]='S' then swds[6]:='-' else swds[6]:='+';
           close(dt);

           assign(dt,'ORBTR2A.sv');
           reset(dt);
           abpq:=-1;
           while not(eof(dt)) do
           begin
                read(dt,swds2,ch,comp10a);
                if (swds2=swds) and (comp10a=comp10) then
                begin
                      read(dt,abpq,evq,thq,prdq,aiq,ommq,ombq,
                         t0q,
                         iro0q,iteta0q,rotq,
                         tetatq,rq,vrq,massaq,f0q,ntq,iq);
                end;
                readln(dt);

          end;
          close(dt);

           assign(dt,'ORBTR2b.sv');
           reset(dt);
           abpw:=-1;
           while not(eof(dt)) do
           begin
                read(dt,swds2,ch,comp10a);
                if (swds2=swds) and (comp10a=comp10) then
                begin
                      read(dt,abpw,evw,thw,prdw,aiw,ommw,ombw,
                         t0w,
                         iro0w,iteta0w,rotw,
                         tetatw,rw,vrw,massaw,f0w,ntw,iw);
                end;
                readln(dt);

          end;
          close(dt);

           assign(dt,'ORBTR1.sv');
           reset(dt);
           abpy:=-1;
           while not(eof(dt)) do
           begin
                read(dt,swds2,ch,comp10a);
                if (swds2=swds) then
                begin
                      compy:=comp10a;
                      read(dt,abpy,evy,thy,prdy,aiy,ommy,omby,
                         t0y,
                         iro0y,iteta0y,roty,
                         tetaty,ry,vry,massay,f0y,nty,iy);
                end;
                readln(dt);

          end;
          close(dt);


          write(rs,swds,' ',comp10);
          if (f0w-f0q)>oev90 then
          begin
               write(rs,'No  ');
               ch:='N';
          end else
          begin
               write(rs,'Yes ');
               ch:='Y';
          end;
          write(rs,f0q:8:5,' ',f0q+oev90:8:5,' ',f0w:8:5,' ',(f0w-f0q)/oev90:7:3);
          writeln(rs);
          ug1:=GetUgol2(abpy,evy,thy,prdy,aiy,omby,ommy,
                        abpq,evq,thq,prdq,aiq,ommq,ombq);
          ug2:=GetUgol2(abpy,evy,thy,prdy,aiy,omby,ommy,
                        abpw,evw,thw,prdw,aiw,ommw,ombw);
          ug1:=ug1*180/pi;
          ug2:=ug2*180/pi;
          if omby>pi then
          begin
               omby:=omby-pi;
               ommy:=ommy-pi;
               if ommy<0 then ommy:=2*pi+ommy;
          end;

          if ombq>pi then
          begin
               ombq:=ombq-pi;
               ommq:=ommq-pi;
               if ommq<0 then ommq:=2*pi+ommq;
          end;

          if ombw>pi then
          begin
               ombw:=ombw-pi;
               ommw:=ommw-pi;
               if ommw<0 then ommw:=2*pi+ommw;
          end;
          wcomp:=copy(compy,1,5);

          writeln(rs2,swds,' ',wcomp:5,' ','     ',' ',f0y:8:5,' ',
          prdy:8:1,' ',abpy:7:3,' ',(pi-aiy)*180/pi:7:2,' ',omby*180/pi:7:2
          ,' ',thy:6:1,' ',evy:6:4,' ',ommy*180/pi:7:2);

          wcomp:=copy(comp10,1,5);
          writeln(rs2,swds,' ',wcomp:5,'  ',ug1:4:1,' ',f0q:8:5,' ',
          prdq:8:1,' ',abpq:7:3,' ',(pi-aiq)*180/pi:7:2,' ',ombq*180/pi:7:2
          ,' ',thq:6:1,' ',evq:6:4,' ',ommq*180/pi:7:2);

          writeln(rs2,swds,' ',wcomp:5,' ',ch,ug2:4:1,' ',f0w:8:5,' ',
          prdw:8:1,' ',abpw:7:3,' ',(pi-aiw)*180/pi:7:2,' ',ombw*180/pi:7:2
          ,' ',thw:6:1,' ',evw:6:4,' ',ommw*180/pi:7:2);
          writeln(rs2);
          i:=0;


      end;
      close(rs);
      close(rs2);

end;


begin
      xfk[ 1]:=-2.20038592302443E+0000;
      xfk[ 2]:= 3.95320923444222E-0001;
      xfk[ 3]:=-1.40651047392788E-0002;
      xfk[ 4]:= 2.97817055140045E-0004;
      xfk[ 5]:=-3.75276309676466E-0006;
      xfk[ 6]:= 2.96148505471758E-0008;
      xfk[ 7]:=-1.51821284077818E-0010;
      xfk[ 8]:= 5.13531862870052E-0013;
      xfk[ 9]:=-1.13748354489786E-0015;
      xfk[10]:= 1.58736519209883E-0018;
      xfk[11]:=-1.26633323012809E-0021;
      xfk[12]:= 4.40435754835436E-0025;

end.