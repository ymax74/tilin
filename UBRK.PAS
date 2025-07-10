{$A+,B-,D+,E-,F+,G+,I+,L+,N+,O+,P+,Q+,R+,S+,T-,V+,X+,Y+}
{$M 65520,0,655360}
unit ubrk ;
interface
uses mnk,mpdq,mpvd,sbop,crt;
type
    ex50=array[1..50] of extended;
    procedure tcir(ti0,tiro0,tiro,titeta0,titeta,
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

    procedure tcir2(ti0,tiro0,tiro,titeta0,titeta,
               tirot,tirot0,titetat0,titetat,
               tir0,tir,tivr0,tivr,tim0,tim,tip0,tip
               :das;
               var rog0,tetag0,eps,raz:extended;
               var i:integer;
               var jb,id:integer;
               et,ro,teta:dst3;
               var rog,tetag,rog1,tetag1,dro:extended;
               bro0,bteta0,brot,btetat,br,bvr,bmassa,bpia:boolean;
               a:msls2);
          procedure ogr1(var raz:extended);
          procedure ogr2(var raz:extended);
          procedure erval(psd:string; code:integer);
          procedure tilin(var l:extended; biskl:boolean;

                nt,t3:integer;
                gbop:idst3;
                ti0:das;
                et,x,y:dst3;
                a:msls);
          procedure ishdata2(pst:string);
          procedure ishdata3(pst:string);
          procedure ish26r;
          procedure ish26r3;
          procedure getves;
          procedure getpvd;
          procedure GetVneIskl;


          var
          t50000:longint;
          chm,chb:string[26];



implementation
procedure  getgaia;
var dt:text;
lwds:string[10];
ch:char;
lcomp5:string[5];
lteta,lro:extended;

begin
     exit;
     assign(dt,'gaia.txt');
     reset(dt);
     rog:=-1;
     while not(eof(dt)) do
     begin
          readln(dt,lwds,ch,lcomp5,lteta,lro);
          if (lwds=swds) and (lcomp5=comp10) then
          begin
               rog:=lro;
               tetag:=lteta;
          end;

     end;
     if (nt<t1330) and (rog>0) then
     begin
          nt:=nt+1;
          et^[nt]:=2015;
          ro^[nt]:=rog;
          teta^[nt]:=tetag;
          x^[nt]:=ro^[nt]*sin(teta^[nt]*pi/180);
          y^[nt]:=ro^[nt]*cos(teta^[nt]*pi/180);

          t0:=t0+2015;
          aref^[nt]:='gaia____'
     end;
     close(dt);
end;

procedure ish26r;
var pstk,fl:string;
dt,rss:text;
ref8:string[6];
lep,lte,lro,ra,dc:extended;
lcomp10:string[5];
sr1:string;
ln,i:integer;
begin
                              ra:=(rac)*15*pi/180;
                              dc:=(dcc)*pi/180;
                              ra:=sin(ra);
                              dc:=1/cos(dc);

     pstk:=swds;
     fl:=copy(pstk,1,5);
     if pstk[6]='-' then fl:=fl+'S'
     else fl:=fl+'N';
     fl:=fl+copy(pstk,7,1);
     fl:=fl+'.';
     fl:=fl+copy(pstk,8,3);
     fl:='..\ovd26re2\S'+fl;
     {$I-}
          Assign(dt,fl);
          reset(dt);
     {$I+}
     sr1:='AB';
     sr1:=sr1+' '+' '+' ';

     if IOResult = 0 then
     begin
          while not(eof(dt)) do
          begin
                readln(dt,ref8,lcomp10,lep,lte,lro);
                ln:=length(lcomp10);
                for i:=1 to ln do if lcomp10[i]='-' then lcomp10[i]:=',';
                if (lcomp10=comp10) or
                   ((comp10=sr1) and (lcomp10='     '))
                then
                begin
                     nt:=nt+1;

                     if nt<=t1330 then
                     begin
                                et^[nt]:=lep;
                                teta^[nt]:=lte;
                                ro^[nt]:=lro;
                                aref^[nt]:=ref8+'  ';
                                if (dcc>-100) then
                                begin
                                     teta^[nt]:=teta^[nt]+0.00556*ra*(2000-et^[nt])*dc;
                                end;
                                x^[nt]:=ro^[nt]*sin(teta^[nt]*pi/180);
                                y^[nt]:=ro^[nt]*cos(teta^[nt]*pi/180);


                     end;
                     t0:=t0+lep;

                end;

          end;
          close(dt);
     end;
     getgaia;
          if nt>1 then
          begin
               t0:=t0/nt;
               nnt:=nt;
               if nnt>t1330 then nnt:=t1330;
               for i:=1 to nnt do
               begin
                     et^[i]:=et^[i]-t0;
               end;
          end;
          t4:=1;
          nnt:=nt;
          if nnt>t1330 then nnt:=t1330;

          t1:=t4;
          pia:=1;
          {assign(rss,'sour.dt');
          rewrite(rss);
          writeln(rss,rac:15:9,' ',dcc:14:8);
          for i:=1 to nnt do
          begin
               writeln(rss,et^[i]+t0:10:5,' ',teta^[i]:10:5,' ',ro^[i]:10:5,' ',aref^[i]);
          end;
          close(rss);}



end;


procedure sort(tep:ex50;kt:integer);
var ii,jj,i:integer;
x:extended;
hy:boolean;
begin
     i:=kt;
     for ii:=2 to i do
     begin
          for jj:=i downto ii do
          begin
               hy:=(tep[jj-1]>tep[jj]);
               if hy then
               begin
                     x:=tep[jj-1];
                     tep[jj-1]:=tep[jj];
                     tep[jj]:=x;
               end;
          end;
      end;
 end;

procedure sort2(tep:dst2;kt:integer);
var ii,jj,i:integer;
x:extended;
hy:boolean;
begin
     i:=kt;
     for ii:=2 to i do
     begin
          for jj:=i downto ii do
          begin
               hy:=(tep^[jj-1]<tep^[jj]);
               if hy then
               begin
                     x:=tep^[jj-1];
                     tep^[jj-1]:=tep^[jj];
                     tep^[jj]:=x;
               end;
          end;
      end;
 end;

function otnct(rf:string):string;
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
     k:=length(pst);
     if pst[k]='.' then pst:=copy(pst,1,k-1);

     otnct:=pst;

end;

procedure ishdata2(pst:string);
var dt:text;
rss:text;
vzv,vzv2:boolean;
i,er1,er2,er3,er4,er5,er6,nsj,tncomp,gu:integer;
lep,lro,lte:extended;
pst2,pst3:string;
ra,rh,rs,dc,dg,dm,ds:extended;


begin
          ecm:=true;
          assign(dt,pst);
          reset(dt);
          vzv:=false;
          vzv2:=false;
          nt:=0;
          dc:=-100;
          rac:=0;
          dcc:=-100;
          t0:=0;
          nsj:=0;
          tncomp:=1;
          while not(eof(dt)) and not(vzv2) do
          begin
               nt:=nt+1;
               readln(dt, lep,lte,lro);
                           if nt<=t1330 then
                           begin
                                et^[nt]:=lep;
                                teta^[nt]:=lte;
                                ro^[nt]:=lro;
                                aref^[nt]:='--------';
                                x^[nt]:=ro^[nt]*sin(teta^[nt]*pi/180);
                                y^[nt]:=ro^[nt]*cos(teta^[nt]*pi/180);
                           end;
                           t0:=t0+lep;


          end;
          exit;
          while not(eof(dt)) and not(vzv2) do
          begin
               nsj:=nsj+1;
               pst:='';
               ch:=' ';
               gu:=0;
               while (ch<>#10) and (gu<256) do
               begin
                    read(dt,ch);
                    pst:=pst+ch;
                    gu:=gu+1;
               end;
               if nsj=1 then
               begin
                    swds:=copy(pst,31,10);
               end;
               i:=(length(pst));

               if vzv
                  and  (pst[7]=' ')
                  and  ('0'<=pst[11]) and (pst[11]<='9')
                  and  (pst[23]='.')
                  and  (pst[39]='.')
                  and  (length(pst)>39)

               then
               begin
                      pst2:=copy(pst,8,11);
                      pst2:=otnct(pst2);
                      val(pst2,lep,er1);

                      pst2:=copy(pst,20,7);
                      pst2:=otnct(pst2);
                      val(pst2,lte,er2);

                      pst2:=copy(pst,36,8);
                      pst2:=otnct(pst2);
                      val(pst2,lro,er3);
                      if (er1=0) and (er2=0) and (er3=0) then
                      begin
                           nt:=nt+1;
                           if nt<=t1330 then
                           begin
                                et^[nt]:=lep;
                                teta^[nt]:=lte;
                                ro^[nt]:=lro;
                                aref^[nt]:=copy(pst,103,8);
                                if (dc>-100) then
                                begin
                                     teta^[nt]:=teta^[nt]+0.00556*ra*(2000-et^[nt])*dc;
                                end;
                                x^[nt]:=ro^[nt]*sin(teta^[nt]*pi/180);
                                y^[nt]:=ro^[nt]*cos(teta^[nt]*pi/180);
                           end;
                           t0:=t0+lep;


                      end;
                      i:=0;

               end else
               begin
                    if vzv then vzv2:=true;
               end;
               if      (pst[1]<>' ')
                  and  ('0'<=pst[27]) and (pst[27]<='9')
                  and  (pst[28]=' ')
                  and  (pst[109]='.')
                  and  (pst[119]='.')
                  and  (length(pst)>119)
                  and  not(vzv)
                  then
               begin
                    if tncomp=ncomp then
                    begin
                         vzv:=true;
                         comp10:=copy(pst,8,5);

                         pst3:=copy(pst,103,18);

                         pst2:=copy(pst3,1,2);
                         val(pst2,rh,er1);

                         pst2:=copy(pst3,3,2);
                         val(pst2,rm,er2);

                         pst2:=copy(pst3,5,5);
                         val(pst2,rs,er3);

                         pst2:=copy(pst3,11,2);
                         val(pst2,dg,er4);

                         pst2:=copy(pst3,13,2);
                         val(pst2,dm,er5);

                         pst2:=copy(pst3,15,4);
                         val(pst2,ds,er6);

                         if (er1=0) and (er2=0) and (er3=0) and
                            (er4=0) and (er5=0) and (er6=0) then
                         begin
                              rac:=(rh+rm/60+rs/3600.0);
                              dcc:=(dg+dm/60+ds/3600);
                              ra:=(rh+rm/60+rs/3600.0)*15*pi/180;
                              dc:=(dg+dm/60+ds/3600)*pi/180;
                              if pst3[10]='-' then
                              begin
                                   dc:=-dc;
                                   dcc:=-dcc;
                              end;
                              ra:=sin(ra);
                              dc:=1/cos(dc);
                         end;
                    end else
                    begin
                         tncomp:=tncomp+1;
                    end;
               end;

          end;
          ecm:=not(eof(dt));
          close(dt);

end;

procedure ish26r3;
var pstk,fl:string;
dt,rss:text;
ref8:string[6];
lep,lte,lro,srdx,srdy,min,ras,ra,dc:extended;
lcomp10:string[5];
sr1:string;
ln,i,tnt,mi,er1,er2:integer;
alep,alro,alte:ex50;
ax,ay:ex50;
alref:array[1..50] of string[8];
vyh:boolean;

begin
                              ra:=(rac)*15*pi/180;
                              dc:=(dcc)*pi/180;
                              ra:=sin(ra);
                              dc:=1/cos(dc);

          i:=mnd;
          er1:=trunc(int(i/2));
          er2:=i-er1*2;
          vyh:=er2=0;

     pstk:=swds;
     fl:=copy(pstk,1,5);
     if pstk[6]='-' then fl:=fl+'S'
     else fl:=fl+'N';
     fl:=fl+copy(pstk,7,1);
     fl:=fl+'.';
     fl:=fl+copy(pstk,8,3);
     fl:='..\ovd26re2\S'+fl;
     {$I-}
          Assign(dt,fl);
          reset(dt);
     {$I+}
     sr1:='AB';
     sr1:=sr1+' '+' '+' ';
     tnt:=0;

     if IOResult = 0 then
     begin
          while not(eof(dt)) do
          begin
                readln(dt,ref8,lcomp10,lep,lte,lro);
                ln:=length(lcomp10);
                for i:=1 to ln do if lcomp10[i]='-' then lcomp10[i]:=',';
                if (lcomp10=comp10) or
                   ((comp10=sr1) and (lcomp10='     '))
                then
                begin
                           tnt:=tnt+1;
                           {if nt<=t1330 then}
                           begin
                                alep[tnt]:=lep;
                                alte[tnt]:=lte;
                                alro[tnt]:=lro;
                                alref[tnt]:=ref8;
                                if (dcc>-99) then
                                begin
                                     alte[tnt]:=alte[tnt]+0.00556*ra*(2000-alep[tnt])*dc;
                                end;
                                ax[tnt]:=alro[tnt]*sin(alte[tnt]*pi/180);
                                ay[tnt]:=alro[tnt]*cos(alte[tnt]*pi/180);
                           end;

                           {if nt<=t1330 then}
                           if tnt=mnd then
                           begin
                                sort(ax,tnt);
                                sort(ay,tnt);
                                if vyh then
                                begin
                                     mi:=trunc(tnt/2);
                                     srdx:=(ax[mi]+ax[mi+1])*0.5;
                                     srdy:=(ay[mi]+ay[mi+1])*0.5;
                                end else
                                begin
                                     mi:=trunc(tnt/2)+1;
                                     srdx:=ax[mi];
                                     srdy:=ay[mi];
                                end;
                                mi:=-1;
                                min:=1e39;
                                for j:=1 to tnt do
                                begin
                                     ras:=sqrt(sqr(ax[j]-srdx)+sqr(ay[j]-srdy));
                                     if ras<min then
                                     begin
                                          min:=ras;
                                          mi:=j;
                                     end;
                                end;
                                nt:=nt+1;
                                et^[nt]:=alep[mi];
                                teta^[nt]:=alte[mi];
                                ro^[nt]:=alro[mi];
                                aref^[nt]:=alref[mi]+'  ';
                                x^[nt]:=ax[mi];
                                y^[nt]:=ay[mi];
                                tnt:=0;
                                t0:=t0+et^[nt];

                           end;

                end;

          end;
          close(dt);
     end;
     getgaia;
          if nt>1 then
          begin
               t0:=t0/nt;
               nnt:=nt;
               if nnt>t1330 then nnt:=t1330;
               for i:=1 to nnt do
               begin
                     et^[i]:=et^[i]-t0;
               end;
          end;
          t4:=1;
          nnt:=nt;
          if nnt>t1330 then nnt:=t1330;

          t1:=t4;
          pia:=1;
          {assign(rss,'sour.dt');
          rewrite(rss);
          writeln(rss,rac:15:9,' ',dcc:14:8);
          for i:=1 to nnt do
          begin
               writeln(rss,et^[i]+t0:10:5,' ',teta^[i]:10:5,' ',ro^[i]:10:5,' ',aref^[i]);
          end;
          close(rss);}



end;

procedure ishdata3(pst:string);
var dt:text;
rss:text;
vzv,vzv2:boolean;
i,er1,er2,er3,er4,er5,er6,nsj,tncomp,tnt,mi,j:integer;
lep,lro,lte:extended;
alep,alro,alte:ex50;
ax,ay:ex50;
pst2,pst3:string;
ra,rh,rs,dc,dg,dm,ds,srx:extended;
vyh,vyh2:boolean;
alref:array[1..50] of string[8];
srdx,srdy,min,ras:extended;



begin
          ecm:=true;
          i:=2;
          {nt:=751;}
          repeat
                i:=i+1;
                er1:=trunc(int(nt/i));
                er2:=nt-er1*i;
                vyh:=(er2=0) and (er1<=t1330);
                vyh2:=(er1+1)<=t1330;
                {er2:=0;}

          until vyh or vyh2;
          {i:=5;}
          mnd:=i;
          er1:=trunc(int(i/2));
          er2:=i-er1*2;
          vyh:=er2=0;

          assign(dt,pst);
          reset(dt);
          vzv:=false;
          vzv2:=false;
          nt:=0;
          dc:=-100;
          rac:=0;
          dcc:=-100;
          t0:=0;
          nsj:=0;
          tncomp:=1;
          tnt:=0;
          while not(eof(dt)) and not(vzv2) do
          begin
               nsj:=nsj+1;
               {pst:='';
               ch:=' ';
               while ch<>#10 do
               begin
                    read(dt,ch);
                    pst:=pst+ch;
               end;
               if nsj=1 then
               begin
                    swds:=copy(pst,31,10);
               end;
               i:=(length(pst));

               if vzv
                  and  (pst[7]=' ')
                  and  ('0'<=pst[11]) and (pst[11]<='9')
                  and  (pst[23]='.')
                  and  (pst[39]='.')
                  and  (length(pst)>39)

               then}
               begin
                      {pst2:=copy(pst,8,11);
                      pst2:=otnct(pst2);
                      val(pst2,lep,er1);

                      pst2:=copy(pst,20,7);
                      pst2:=otnct(pst2);
                      val(pst2,lte,er2);

                      pst2:=copy(pst,36,8);
                      pst2:=otnct(pst2);
                      val(pst2,lro,er3);
                      if (er1=0) and (er2=0) and (er3=0) then}
                      readln(dt,lep,lte,lro);
                      begin

                           tnt:=tnt+1;
                           {if nt<=t1330 then}
                           begin
                                alep[tnt]:=lep;
                                alte[tnt]:=lte;
                                alro[tnt]:=lro;
                                alref[tnt]:='--------';
                                {if (dc>-100) then
                                begin
                                     alte[tnt]:=alte[tnt]+0.00556*ra*(2000-alep[tnt])*dc;
                                end;}
                                ax[tnt]:=alro[tnt]*sin(alte[tnt]*pi/180);
                                ay[tnt]:=alro[tnt]*cos(alte[tnt]*pi/180);
                           end;

                           {if nt<=t1330 then}
                           if tnt=mnd then
                           begin
                                sort(ax,tnt);
                                sort(ay,tnt);
                                if vyh then
                                begin
                                     mi:=trunc(tnt/2);
                                     srdx:=(ax[mi]+ax[mi+1])*0.5;
                                     srdy:=(ay[mi]+ay[mi+1])*0.5;
                                end else
                                begin
                                     mi:=trunc(tnt/2)+1;
                                     srdx:=ax[mi];
                                     srdy:=ay[mi];
                                end;
                                mi:=-1;
                                min:=1e39;
                                for j:=1 to tnt do
                                begin
                                     ras:=sqrt(sqr(ax[j]-srdx)+sqr(ay[j]-srdy));
                                     if ras<min then
                                     begin
                                          min:=ras;
                                          mi:=j;
                                     end;
                                end;
                                nt:=nt+1;
                                et^[nt]:=alep[mi];
                                teta^[nt]:=alte[mi];
                                ro^[nt]:=alro[mi];
                                aref^[nt]:=alref[mi];
                                x^[nt]:=ax[mi];
                                y^[nt]:=ay[mi];
                                tnt:=0;
                                t0:=t0+et^[nt];

                           end;


                      end;
                      i:=0;

               end;

               {else
               begin
                    if vzv then vzv2:=true;
               end;
               if      (pst[1]<>' ')
                  and  ('0'<=pst[27]) and (pst[27]<='9')
                  and  (pst[28]=' ')
                  and  (pst[109]='.')
                  and  (pst[119]='.')
                  and  (length(pst)>119)
                  and  not(vzv)
                  then
               begin
                    if tncomp=ncomp then
                    begin
                         vzv:=true;
                         comp10:=copy(pst,8,5);

                         pst3:=copy(pst,103,18);

                         pst2:=copy(pst3,1,2);
                         val(pst2,rh,er1);

                         pst2:=copy(pst3,3,2);
                         val(pst2,rm,er2);

                         pst2:=copy(pst3,5,5);
                         val(pst2,rs,er3);

                         pst2:=copy(pst3,11,2);
                         val(pst2,dg,er4);

                         pst2:=copy(pst3,13,2);
                         val(pst2,dm,er5);

                         pst2:=copy(pst3,15,4);
                         val(pst2,ds,er6);

                         if (er1=0) and (er2=0) and (er3=0) and
                            (er4=0) and (er5=0) and (er6=0) then
                         begin
                              rac:=(rh+rm/60+rs/3600.0);
                              dcc:=(dg+dm/60+ds/3600);
                              ra:=(rh+rm/60+rs/3600.0)*15*pi/180;
                              dc:=(dg+dm/60+ds/3600)*pi/180;
                              if pst3[10]='-' then
                              begin
                                   dc:=-dc;
                                   dcc:=-dcc;
                              end;
                              ra:=sin(ra);
                              dc:=1/cos(dc);
                         end;
                    end else
                    begin
                         tncomp:=tncomp+1;
                    end;
               end;}

          end;
          ecm:=not(eof(dt));
          close(dt);

                           {tnt:=5;}
                           if tnt>0 then
                           begin
                                sort(ax,tnt);
                                sort(ay,tnt);
                                i:=tnt mod 2;
                                if i=0 then
                                begin
                                     mi:=trunc(tnt/2);
                                     srdx:=(ax[mi]+ax[mi+1])*0.5;
                                     srdy:=(ay[mi]+ay[mi+1])*0.5;
                                end else
                                begin
                                     mi:=trunc(tnt/2)+1;
                                     srdx:=ax[mi];
                                     srdy:=ay[mi];
                                end;
                                mi:=-1;
                                min:=1e39;
                                for j:=1 to tnt do
                                begin
                                     ras:=sqrt(sqr(ax[j]-srdx)+sqr(ay[j]-srdy));
                                     if ras<min then
                                     begin
                                          min:=ras;
                                          mi:=j;
                                     end;
                                end;
                                nt:=nt+1;
                                et^[nt]:=alep[mi];
                                teta^[nt]:=alte[mi];
                                ro^[nt]:=alro[mi];
                                aref^[nt]:=alref[mi];
                                x^[nt]:=ax[mi];
                                y^[nt]:=ay[mi];
                                tnt:=0;
                                t0:=t0+et^[nt];
                           end;


end;
procedure ogr1(var raz:extended);
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
end;

procedure tcir(ti0,tiro0,tiro,titeta0,titeta,
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
                     jb:=0;
                     {i:=id;}
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
                     {a^[i*2,jb+1]:=raz*pi*rog/180;}
end;

procedure tcir2(ti0,tiro0,tiro,titeta0,titeta,
               tirot,tirot0,titetat0,titetat,
               tir0,tir,tivr0,tivr,tim0,tim,tip0,tip
               :das;
               var rog0,tetag0,eps,raz:extended;
               var i:integer;
               var jb,id:integer;
               et,ro,teta:dst3;
               var rog,tetag,rog1,tetag1,dro:extended;
               bro0,bteta0,brot,btetat,br,bvr,bmassa,bpia:boolean;
               a:msls2);
begin
                     jb:=0;
                     {i:=id;}
                     i:=i+1;
                     mpdq.rtti(ti0,et^[id],rog,tetag);
                     if bro0 then
                     begin
                          mpdq.rtti(tiro0,et^[id],rog0,tetag0);
                          mpdq.rtti(tiro,et^[id],rog1,tetag1);
                          a^[i*2-1]^[jb+1]:=aves^[id]*(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=aves^[id]*raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;

                     if bteta0 then
                     begin
                          mpdq.rtti(titeta0,et^[id],rog0,tetag0);
                          mpdq.rtti(titeta,et^[id],rog1,tetag1);
                          dro:=eps*180/pi;
                          a^[i*2-1]^[jb+1]:=aves^[id]*(rog1-rog0)/dro;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=aves^[id]*raz*pi*ro^[id]/(dro*180);
                          jb:=jb+1;
                     end;
                     if brot then
                     begin
                          mpdq.rtti(tirot0,et^[id],rog0,tetag0);
                          mpdq.rtti(tirot,et^[id],rog1,tetag1);
                          dro:=eps*0.01;
                          a^[i*2-1]^[jb+1]:=aves^[id]*(rog1-rog0)/dro;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=aves^[id]*raz*pi*ro^[id]/(dro*180);
                          jb:=jb+1;
                     end;
                     if btetat then
                     begin
                          mpdq.rtti(titetat0,et^[id],rog0,tetag0);
                          mpdq.rtti(titetat,et^[id],rog1,tetag1);
                          a^[i*2-1]^[jb+1]:=aves^[id]*(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=aves^[id]*raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;
                     if br then
                     begin
                          mpdq.rtti(tir0,et^[id],rog0,tetag0);
                          mpdq.rtti(tir,et^[id],rog1,tetag1);
                          a^[i*2-1]^[jb+1]:=aves^[id]*(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=aves^[id]*raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;


                     if bvr then
                     begin
                          mpdq.rtti(tivr0,et^[id],rog0,tetag0);
                          mpdq.rtti(tivr,et^[id],rog1,tetag1);
                          a^[i*2-1]^[jb+1]:=aves^[id]*(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=aves^[id]*raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;
                     if bmassa then
                     begin
                          mpdq.rtti(tim0,et^[id],rog0,tetag0);
                          mpdq.rtti(tim,et^[id],rog1,tetag1);
                          a^[i*2-1]^[jb+1]:=aves^[id]*(rog1-rog0)/eps;
                          raz:=(tetag1-tetag0);
                          ogr1(raz);
                          a^[i*2]^[jb+1]:=aves^[id]*raz*pi*ro^[id]/(eps*180);
                          jb:=jb+1;
                     end;

                     if bpia then
                     begin
                          mpdq.rtti(tip0,et^[id],rog0,tetag0);
                          mpdq.rtti(tip,et^[id],rog1,tetag1);
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
                     {a^[i*2,jb+1]:=raz*pi*rog/180;}
end;
function fnev2(nt:integer;ti0:das;gbop:idst3;et,x,y:dst3):extended;
var
dro,dte,raz,sum,rog,tetag,xg,yg,dx,dy,eb,ceb,seb,xx,yy:extended;
nt2,id:integer;
begin

      sum:=0;
      nt2:=0;
      for id:=1 to nt do
      begin
           if gbop^[id]=0 then
           begin
                nt2:=nt2+1;
                mpdq.rtti2(ti0,et^[id],xg,yg,eb,ceb,seb,xx,yy);
                dx:=xg-x^[id];
                dy:=yg-y^[id];
                sum:=sum+dx*dx+dy*dy;
           end;
      end;
      sum:=0.5*sum/nt2;
      sum:=sqrt(sum);
      fnev2:= sum;
end;
function flam2(lam:extended;vgrad:das;nt:integer;ti0:das;gbop:idst3;et,x,y:dst3):extended;
var    ti:das;

begin
                     {jb:=0;
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
                     end;}
                     ti[1]:=ti0[1];
                     ti[2]:=ti0[2]-lam*vgrad[1];
                     ti[3]:=ti0[3]-lam*vgrad[2];
                     ti[4]:=ti0[4]-lam*vgrad[3];
                     ti[5]:=ti0[5]-lam*vgrad[4];
                     ti[6]:=ti0[6]-lam*vgrad[5];
                     ti[7]:=ti0[7]-lam*vgrad[6];
                     ti[8]:=ti0[8]-lam*vgrad[7];



        flam2:=fnev2(nt,ti,gbop,et,x,y);

end;


procedure tilin(var l:extended; biskl:boolean;
                nt,t3:integer;
                gbop:idst3;
                ti0:das;
                et,x,y:dst3;
                a:msls);

var ih,id2:integer;
    lp,rog2,ad,bd,h:extended;
    rsx:text;
    vgrad:das;

    ch:char;
    f0,f1,f2,f3,f4,maxnev,ff1,ff2:extended;
    mi,itr:integer;
    sbvr,sbmassa,sbpia:boolean;
    spia,svmassa,svvr,ln:extended;
    uge:boolean;
    siro0,siteta0,srot,stetat,sr,eps2,eps,svll,svlll,ev,xq,yq,eb,ceb,seb,ebt,ebe,ebp
    ,mm,tt,pp,aa,bb,ff,gg,xx,yy,s1ev,q1,q2,xq2,yq2,pg,eb2,ro,tetag,teta,
    brt,btt,dx,dy,rof,tetaf,raz,ro2,teta2:extended;

    rsh:text;
    fd:longint;
    xf,sgm:stb;
    vyh:boolean;
    ti,ti2:das;
    i:integer;
    jb,tgu:integer;





begin
      {ti0[5]:=0.1;}
      {cls;}
      ClrScr;
      eps:=1e-7;
      eps2:=1e-6;
      l:=1e-5;

      {bvr:=false;
      bmassa:=false;
      bpia:=false;
      pia:=1;}
      {br:=false;
      bro0:=false;
      bteta0:=false;}
      {brot:=false;
      btetat:=false;}


      ch:='p';

          fd:=0;
          for ih:=1 to 8 do ti[ih]:=ti0[ih];
          ti[5]:=ti[1]*ti[5];
          ti[6]:=ti[1]*ti[6];
          ti[7]:=ti[1]*ti[7];
          ti[8]:=ti[1]*ti[8];
          ti[1]:=1;


          begin
                {massa:=massa-0.1;}
          repeat
          t3:=nt;


          repeat
                ev:=ti0[2];
                i:=0;
                jb:=0;
                for id2:=1 to t3 do
                begin
                     if gbop^[id2]=0 then
                     begin
                          jb:=jb+1;
                     end;
                end;
                for id2:=1 to t3 do
                begin
                     if gbop^[id2]=0 then
                     begin
                          i:=i+1;
                          mpdq.rtti2(ti,et^[id2],xq,yq,eb,ceb,seb,xx,yy);
                          {ebt:=(et^[id2]-1)*2*pi/((1-ti[2]*ceb)*ti[4]);}
                          ebt:=(-2*pi/(ti[4]))/(1-ti[2]*ceb);
                          mm:=(et^[id2]-ti[3])*2*pi/ti[4];
                          ebe:=(seb)/(1-ti[2]*ceb);
                          ebp:=-mm/(ti[4]*(1-ti[2]*ceb));
                          ev:=ti[2];
                          tt:=ti[3];
                          pp:=ti[4];
                          aa:=ti[5];
                          bb:=ti[6];
                          ff:=ti[7];
                          gg:=ti[8];
                          s1ev:=sqrt(1-ev*ev);
                          q1:=ff*s1ev*ceb-aa*seb;
                          q2:=-bb*seb+gg*s1ev*ceb;
                          ro:=sqrt(xq*xq+yq*yq);
                          tetag:=arctan(xq/yq);
                          if yq<0 then tetag:=pi+tetag;
                          if tetag<0 then tetag:=tetag+2*pi;
                          teta:=tetag*180/pi;

                          rof:=sqrt(sqr(x^[id2])+sqr(y^[id2]));

                          tetaf:=arctan(x^[id2]/y^[id2]);
                          if y^[id2]<0 then tetaf:=pi+tetaf;
                          if tetaf<0 then tetaf:=tetaf+2*pi;
                          tetaf:=tetaf*180/pi;


                          {for ih:=1 to 8 do ti2[ih]:=ti[ih];
                          tgu:=7;
                          ti2[tgu+1]:=ti2[tgu+1]+eps;
                          mpdq.rtti2(ti2,et^[id2],xq2,yq2,eb2,ceb,seb,xx,yy);
                          pg:=(yq2-yq)/eps;
                          ro2:=sqrt(xq2*xq2+yq2*yq2);
                          teta2:=arctan(xq2/yq2);
                          if yq2<0 then teta2:=pi+teta2;
                          if teta2<0 then teta2:=teta2+2*pi;
                          teta2:=teta2*180/pi;
                          pg:=rof*(teta2-teta)*pi/(eps*180);}




                          {brt:=1/(2*ro);
                          btt:=rof/(1+sqr(xq/yq));

                          dx:=q1*ebe-(aa+ff*ev*seb/s1ev);
                          dy:=q2*ebe-(bb+gg*ev*seb/s1ev);
                          a^[i*2-1,1]:=brt*(2*xq*dx+2*yq*dy);
                          a^[i*2  ,1]:=btt*((dx/yq)-dy*xq/(yq*yq));

                          dx:=q1*ebt;;
                          dy:=q2*ebt;;
                          a^[i*2-1,2]:=brt*(2*xq*dx+2*yq*dy);
                          a^[i*2  ,2]:=btt*((dx/yq)-dy*xq/(yq*yq));

                          dx:=q1*ebp;;
                          dy:=q2*ebp;;
                          a^[i*2-1,3]:=brt*(2*xq*dx+2*yq*dy);
                          a^[i*2  ,3]:=btt*((dx/yq)-dy*xq/(yq*yq));

                          dx:=xx;
                          dy:=0;
                          a^[i*2-1,4]:=brt*(2*xq*dx+2*yq*dy);
                          a^[i*2  ,4]:=btt*((dx/yq)-dy*xq/(yq*yq));

                          dx:=0;
                          dy:=xx;
                          a^[i*2-1,5]:=brt*(2*xq*dx+2*yq*dy);
                          a^[i*2  ,5]:=btt*((dx/yq)-dy*xq/(yq*yq));

                          dx:=yy;
                          dy:=0;
                          a^[i*2-1,6]:=brt*(2*xq*dx+2*yq*dy);
                          a^[i*2  ,6]:=btt*((dx/yq)-dy*xq/(yq*yq));

                          dx:=0;
                          dy:=yy;
                          a^[i*2-1,7]:=brt*(2*xq*dx+2*yq*dy);
                          a^[i*2  ,7]:=btt*((dx/yq)-dy*xq/(yq*yq));

                          a^[i*2-1,8]:=ro-rof;
                          raz:=(teta-tetaf);
                          ogr2(raz);
                          a^[i*2  ,8]:=raz*pi*rof/180;;}


                          a^[i*2-1,1]:=q1*ebe-(aa+ff*ev*seb/s1ev);
                          a^[i*2-1,2]:=q1*ebt;
                          a^[i*2-1,3]:=q1*ebp;
                          a^[i*2-1,4]:=xx;
                          a^[i*2-1,5]:=0;
                          a^[i*2-1,6]:=yy;
                          a^[i*2-1,7]:=0;

                          a^[i*2  ,1]:=q2*ebe-(bb+gg*ev*seb/s1ev);
                          a^[i*2  ,2]:=q2*ebt;
                          a^[i*2  ,3]:=q2*ebp;
                          a^[i*2  ,4]:=0;
                          a^[i*2  ,5]:=xx;
                          a^[i*2  ,6]:=0;
                          a^[i*2  ,7]:=yy;

                          a^[i*2-1,8]:=xq-x^[id2];
                          a^[i*2  ,8]:=yq-y^[id2];


                     end;
                end;
                {halt(0);}
                {close(rsx);}
                ptlst(a,xf,cv,i*2,7);
                {sgmnk(a,cv,xf,sgm,7,i*2);}
                {xf[7]:=0;}

                for ih:=1 to 7 do vgrad[ih]:=xf[ih];


                {f0:=fnev2(nt,ti,gbop,et,x,y);}
                f0:=flam2(0,vgrad,nt,ti,gbop,et,x,y);

                f2:=flam2(l,vgrad,nt,ti,gbop,et,x,y);

                l:=abs(l);
                svll:=l;
                if fd>5000 then
                begin
                while (f2>f0) {and (fd<>258)} do
                begin
                     l:=l/1.1;

                     f2:=flam2(l,vgrad,nt,ti0,gbop,et,x,y);;
                end;
                begin
                     f3:=flam2(l/3,vgrad,nt,ti0,gbop,et,x,y);
                     f4:=flam2(l*2,vgrad,nt,ti0,gbop,et,x,y);
                     if f3<f2 then l:=l/3;
                     if f4<=f2 then l:=l*2;
                end;
                end;


                     jb:=0;
                     vyh:=true;

                     lp:=l;
                     {if bro0 then
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
                     end;}
                     for ih:=2 to 8 do
                     begin
                          ti[ih]:=ti[ih]-lp*xf[ih-1];
                          vyh:=vyh and ((xf[ih-1])<eps);
                     end;
                     fd:=fd+1;
                if fd mod 1 =0 then
                begin
                gotoxy(1,1);
                writeln;
                writeln('     e=',ti[2]:10:5);
                writeln('     T=',ti[3]:10:5);
                writeln('     P=',ti[4]:10:5);
                writeln('    AA=',ti[5]:10:5);
                writeln('    BB=',ti[6]:10:5);
                writeln('    FF=',ti[7]:10:5);
                writeln('    GG=',ti[8]:10:5);


                writeln('iteracija=',fd,' kol-vo urav=',i);

                writeln('f0(oev)=',f0:10:5 );
                writeln('lp=',lp );
                end;
             ch:='p';

            if keypressed then ch:=readkey;

            if (l<1e-10) or (ev>50)then
            begin
                 vyh:=true;
                 l:=1e-15;
            end;
            if ch='a' then halt(0);



          until (ch='a') or vyh  or (fd>t50000);
          if {(l>1e-10) and} (ev<50) then
          begin

          if biskl then
          begin
               {maxiskl(iro0,iteta0,rot,tetat,r,vr,massa,pia,oevr,oevt,maxnev,mi);
               ew:=false;
               if (mi>1) and (mi<nt) then
               begin
                    ew:=maxnev>porog;
                    if ew then
                    begin
                         cvb2(mi);
                         if l<1e-10 then
                         begin
                               iro0:=siro0;
                               iteta0:=siteta0;
                               rot:=srot;
                               tetat:=stetat;
                               r:=sr;
                         end;
                         l:=ln;

                    end;
               end;}
          end;
          end;
          until  (ch='a'){ or not(ew)};
          {if biskl then vydres2(iro0,iteta0,rot,tetat,r,vr,massa,pia,t0,et,ro,teta,nt,oevr,oevt);}
          end;

end;
function toup(pst:string):string;
var i,j:integer;

begin
     for i:=1 to length(pst) do
     begin
          for j:=1 to 26 do
          begin
               if pst[i]=chm[j] then pst[i]:=chb[j];
          end;
     end;
     toup:=pst;
end;
procedure normves;
var i,j,k:integer;
sumves,kvr:extended;

begin
     if nt<6 then exit;
     sumves:=0;
     for i:=1 to nt do
     begin
          sumves:=sumves+aves^[i];
     end;
     if sumves<=0 then exit;
     kvr:=nt/sumves;
     sumves:=0;
     for i:=1 to nt do
     begin
          aves^[i]:=kvr*aves^[i];
          sumves:=sumves+aves^[i];
     end;
     i:=0;



end;

procedure normves2;
var i,j,k:integer;
sumves,kvr,nt2:extended;

begin
     if nt<6 then exit;
     sumves:=0;
     nt2:=0;
     for i:=1 to nt do
     begin
          if gbop^[i]=0 then
          begin
               nt2:=nt2+1;
               sumves:=sumves+aves^[i];
          end;
     end;
     if sumves<=0 then exit;
     kvr:=nt2/sumves;
     sumves:=0;
     for i:=1 to nt do
     begin
          aves^[i]:=kvr*aves^[i];
          sumves:=sumves+aves^[i];
     end;
     i:=0;



end;
procedure ogrves;
var i,j,k:integer;
prves:extended;
begin
      if nt<10 then
      begin
           for i:=1 to nt do
           begin
                aves^[i]:=1;
           end;
           exit;
      end;
      for i:=1 to nt do
      begin
           aex^[i]:=aves^[i];
      end;
      sort2(aex,nt);
      prves:=aex^[10];
      for i:=1 to nt do
      begin
           if aves^[i]>prves then
           begin
                aves^[i]:=prves;
           end;
      end;

      i:=0;
end;

procedure getves;
var dt:text;
i,j,k,mi:integer;
alma:array[1..10] of extended;
pst6:string[6];
ref8:string[8];
ves,oev,v1,ves2:extended;
begin
     if not(bves) then
     begin
          for i:=1 to nt do
          begin
               aves^[i]:=1;
          end;
          exit;

     end;
     for i:=1 to nt do
     begin
          aves^[i]:=-1;
          for j:=1 to 8 do
          begin
               if aref^[i,j]=' ' then aref^[i,j]:='_';
          end;
          aref^[i]:=toup(aref^[i]);
          j:=0;
     end;
     i:=0;
     assign(dt,'vesa3.res');
     {$I-}
     reset(dt);
     {$I+}
     if ioresult<>0 then
     begin
          for i:=1 to nt do
          begin
               aves^[i]:=1;
          end;

          exit;
     end;

     while not(eof(dt)) do
     begin
          read(dt,pst6);
          if pst6='allmal' then
          begin
               read(dt,i,j);
               read(dt,v1,alma[i]);
               alma[i]:=1/(alma[i]*alma[i]);

          end;
          readln(dt);
     end;
     close(dt);

     reset(dt);
     while not(eof(dt)) do
     begin
          readln(dt,ref8,k,v1,ves);
          ref8:=toup(ref8);
          for i:=1 to nt do
          begin
               if ref8=aref^[i] then
               begin
                    if ref8='KIY2001_' then
                    begin
                         v1:=0;
                    end;
                    {writeln(ref8,' ',ves:12:7);}
                    ves2:=1/(ves*ves);
                    aves^[i]:=ves2;
               end;
          end;
          i:=0;
     end;
     close(dt);

     assign(dt,'ref10.res');
     reset(dt);
     while not(eof(dt)) do
     begin
          readln(dt,ref8,k);
          ref8:=toup(ref8);

          for i:=1 to nt do
          begin
               if ref8=aref^[i] then
               begin
                    aves^[i]:=alma[k];
               end;
          end;
          i:=0;
     end;
     for i:=1 to nt do
     begin
               if aves^[i]<0 then
               begin
                    aves^[i]:=alma[1];
               end;
    end;
    close(dt);
    ogrves;
    normves;
    {normves2;}
end;
procedure getpvd;
var
dt:text;
ads5:string[5];
ch:char;
wds10:string[10];
comp2:string[2];
pub,dtvr,svr,mas,pit,lvr:extended;
lcomp:string;
begin
     assign(dt,'vrm3.txt');
     reset(dt);
     readln(dt);
     opit:=0;
     while not(eof(dt)) do
     begin
          readln(dt,ads5,ch,wds10,pub,ch,comp2,dtvr,lvr,svr,mas,pit);
          if comp2='  ' then comp2:='AB';
          lcomp:=comp2;
          while length(lcomp)<5 do lcomp:=lcomp+' ';
          if (swds=wds10) and (lcomp=comp10) then
          begin
               opub:=pub;
               odtvr:=dtvr;
               ovr:=lvr;
               osvr:=svr;
               omas:=mas;
               opit:=pit;
          end;
          ch:=' ';
     end;
     if opit>0.00001 then
     begin
          omas2:=omas;
          ovr2:=ovr;
          omas:=omas*opit*opit*opit/(1000.0*1000.0*1000.0);
          ovr:=ovr/4.74;
          ovr:=ovr*opit/1000;
     end;
     close(dt);

end;
procedure ReadTelo(var dt:text);
var
pst15:string[15];
pst6:string[6];
pst4:string[4];
i,tcomp,mi,j:integer;
ch,ch2:char;
pst8:string[8];
fl2:string;
t02,iro0,iteta0,rot,tetat,r,vr,massa,pia,
glabp,glev,glth,glprd,glomm,glomb,glai,glf0q,
et2,ro2,teta2:extended;




begin
     readln(dt,pst15,t0);
     readln(dt,pst15,iro0);
     readln(dt,pst15,iteta0);
     readln(dt,pst15,rot);
     readln(dt,pst15,tetat);
     readln(dt,pst15,r);
     readln(dt,pst15,vr);
     readln(dt,pst15,massa);
     readln(dt);
     pia:=1;
     readln(dt,pst6,glabp);
     readln(dt,pst6,glev);
     readln(dt,pst6,glth);
     readln(dt,pst6,glprd);
     readln(dt,pst6,glomm);
     readln(dt,pst6,glomb);
     readln(dt,pst6,glai);
     readln(dt);
     readln(dt,pst4,glf0q);
     glai:=180-glai;
     glai:=glai*pi/180;
     glomm:=glomm*pi/180;
     glomb:=glomb*pi/180;
     glth:=glth-t0;

     for i:=1 to 5 do
     begin
          readln(dt);
     end;
     for i:=1 to nt do gbop^[i]:=1;
     repeat
           read(dt,ch);
           if (ch='1') or (ch='2') then
           begin
                 read(dt,et2,ro2,teta2);
                 if ch='1' then et2:=et2+1000;
                 if ch='2' then et2:=et2+2000;
                 et2:=et2-t0;
                 for i:=1 to nt do
                 begin
                      if (abs(et2-et^[i])<1e-5) and (abs(ro2-ro^[i])<1e-5) then
                      begin
                           gbop^[i]:=0;
                      end;

                 end;
           end;
           readln(dt);
     until (ch<>'1') and (ch<>'2');
     for i:=1 to nt do gbop3^[i]:=gbop^[i];
     bVneIskl2:=true;



end;
procedure GetVneIskl;
var dt:text;
i,j,k,tcomp,bu,bu2:integer;
wds10,swds2:string[10];
comp2:string[5];

begin
     spfl('itr001\*.orb',sf,k);
     i:=0;
     bVneIskl2:=false;
     {swds:='1';}
     bu:=0;
     bu2:=0;
     repeat
           i:=i+1;
           assign(dt,'itr001\'+sf^[i]);
           reset(dt);
           repeat
                 read(dt,wds10);
                 if length(wds10)<=7 then readln(dt);
           until (eof(dt)) or (length(wds10)>7);

           if length(wds10)<7 then
           begin
                bu:=bu+1;
           end else
           begin
                if (wds10[6]='+') or (wds10[6]='-') then bu2:=bu2+1;
           end;
           if swds=wds10 then
           begin
           end else
           begin
                close(dt);
           end;
     until (swds=wds10) or (i>=k);
     if swds=wds10 then
     begin
          tcomp:=0;
          j:=0;
          repeat
                if j>0 then read(dt,swds2) else swds2:=wds10;
                readln(dt,ch,comp2);
                j:=j+1;
                if (swds2[6]='+') or (swds2[6]='-') then tcomp:=tcomp+1;
                if (eof(dt)) or (tcomp>15) then
                begin
                     close(dt);
                     exit;
                end;
          until (comp2=comp10) and ((swds2[6]='+') or (swds2[6]='-'));
          if (comp2=comp10) then ReadTelo(dt);
          close(dt);

     end;


end;

begin
     t50000:=50000;
     {assign(dt,'tab.TXT');
     reset(dt);
     for i:=1  to 26 do read(dt,chm[i]);
     readln(dt);
     for i:=1  to 26 do read(dt,chb[i]);
     close(dt);}
     chm:='qwertyuiopasdfghjklzxcvbnm';
     chb:='QWERTYUIOPASDFGHJKLZXCVBNM';


end.
