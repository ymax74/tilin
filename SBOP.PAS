unit sbop;
interface
uses mnk,dos;
const
     {bitrves=true;}
     bitrves=false;
     bves=true;
     {bves=false;}
     {btesorb=true;}
     btesorb=false;
     {bitgtr=true;}
     bitgtr=false;

     boshi=true;
     {boshi=false;}
     bvne=true;
     {bvne=false;}
     bvnevyd=false;
     predmin=1440.0;
               t1330=499;
               {t1330=200;}
               {t1330=26;}
               {t1330=2000;}
               {t16383=16383;}
               t16383=2000;
               t5000=600;
               t2050=2050;
               intver=(1-2*0.34134)/2;
               trs=50;

       (*oprza=true{false};*)
       t100=500;
      {t1330=1580;}

      btorb=false;
      smassa=0.05;
      spia=0.00386;
      svr=0.3;
      porog=5;
      {ovds='..\anorb7';}
      {ovds='..\ovdstst';}
      {ovds='..\ovds05';}
      {ovds='..\ovds';}
      ovds='..\ovdstr2';

      brdln=false;
      t13=15;
      t11=t13-2;
      t8=t13-5;
      {t50000=5000;}
      tt10=20;
      t1024=1;
type pkt=array[1..t5000] of string[12];
     tspfl=^pkt;
     tref=array[1..t1330] of string[6];
     tref2=array[1..t1330] of string[8];
         vector=array[1..3] of extended;
         tdst2=array[1..t1330] of extended;
         dst2=^tdst2;
         itdst2=array[1..t1330] of integer;
         idst2=^itdst2;
         {type} dst=array[1..t1330] of extended;
              idst=array[1..t1330] of integer;
              epk=array[2000..t2050] of extended;
              das=array[1..18] of extended;
              dst3=^dst;
              idst3=^idst;
              ext6=array[1..9] of extended;
              ttrs=array[1..trs] of extended;
              ats2=array[1..t16383] of single;
              ats=^ats2;
        mtreal=extended;






var dt,rpr,rsq,ris,rsv,rso,rsfp:text;
    gbop,gbop2,gbop3,agmi:idst3;
    ngmi:integer;
    ain:idst2;
    x,y,et,ro,aex3,teta{,teta3},aves:dst3;
    {teta:^dst;}
    {x3,y3,et3,ro3:^dst;}
    {arot,atetat,sarot,satetat:^dst;}
    aref{,aref2,aref3}:^tref2;
    uge:array[1..700] of boolean;
    nnt,nnt2,i,ht,sdr,jb,t3,t4,id,nt,t1,ercode,iz,gi,nrd:integer;
    fed,pst:string;
    fd:longint;
    x0,y0,xs,ys,r,om0,omt,omtt,om,cso,sno,l,dxw,dyw,teta1,dfte,eti,t0
    ,teta0,mu,mut,rol,psi,xa,ya,za,xt,yt,zt,kp,rt,abp,ev,ai,omm,omb
    ,fr,sro,dro,dteta,ves,th,prd,kpx,sum,mnro,mnxc,mnyc,gllam:extended;
    comm,somm,comb,somb,si,ci,aa,ff,bb,gg,eb,seb,ceb,xx,yy,xc,yc,rog,tetag
    ,vx,vy,h,ro0,rok,mns,xcp,ycp,omk,ra,dc,svx,svy,e0,via,pc,stt,ctt,tt0,tetat,rot,iteta0,iro0
    ,tetag1,rog1,massa,pia,vr,nro0,sro0,nteta0,steta0,nrot,srot,ntetat,stetat,oevr,oevt,osb,rot2,tetat2,kriv0,skri
    ,sros:extended;
    a,at:msls;
    al:msls2;
    rh,rm,dg,ds,massa0,pi0,le,m1,m2:extended;
    cv:matrls;
    xf,sgm:stb;
    xft,sgmt:stb;
    k,ic:longint;
    bpia,bro0,bteta0,btetat,brot,br,bvr,bmassa,bpi,vyh:boolean;
    ch:char;
    ew,bos,bprr,bnkri,bnvr:boolean;
    aiskl:array[1..300] of integer;
    oevr0,oevt0:extended;
    dtp:text;
    tiv:das;
    mbcir:boolean;
    sf:tspfl;
    ka,k2a,ia,ja:integer;
    comp10,tcomp10:string[10];
    swds:string[10];
    itrf,iud,inud,kf,pkf,sorb1,sorb2,sorb3,sorb4,sorb5,ku,ku2,p1,gmi:integer;
    bprd:boolean;
    tfl,pst2:string;
    aex,aex2:dst2;
    isuge,uci,r26,uci2:boolean;
    gf0,gprd,gnak,tekr:extended;
    ae9:ext6;
    posf:boolean;
    svsgm:array[1..6] of extended;
    mxk2,mdf,dk1,dk2:integer;
    {niskl:integer;}
              {roo,tetao:dst2;}
              fg:boolean;
              brg:boolean;
              xro0,xteta0{roi,tetai,ephi}:^dst;

              otl,rzs,otl2:text;
              abp0,ev0,ai0,omm0,omb0,th0,sabp0,sev0,sai0,somm0,somb0,sth0,prd0,sprd0:extended;
              roe,tetae,sroe,statae:^epk;
              jh,niskl,iv:integer;
              bcir:boolean;
              arsv:array[1..3] of text;
              nevgmi:integer;
              btvht,bzc:boolean;
              mflnm,oflnm:string;
              bmjd:extended;
              rs0:text;
              nevyh,koi:boolean;
              oprza:boolean;
              uabp,uev,uth,uprd,uai,uomm,uomb:boolean;
              pst12:string[12];

rstr,rstr2:text;
awdst:array[1..15] of string[10];
acompt:array[1..15,1..2] of string[5];
atswds:array[1..15] of string[10];
atcomp10:array[1..15] of string[10];
atabp,atev,atth,atprd,atai,atomm,atomb,att0,atiro0,atiteta0,atrot,attetat,
atr,atvr,atmassa,atgf0,atnt,ati:array[1..15] of extended;
fixoi:boolean;
iper:integer;
glev,glth,glprd,glabp,glai,glomb,glomm,glf0q,rac,dcc,gte,gro:extended;
model,fixoi2:boolean;
ncomp,mnd:integer;
ecm:boolean;
opub,odtvr,ovr,osvr,omas,opit,omas2,ovr2:extended;
bklpvd,bVneIskl,bVneIskl2:boolean;







              function tekmjdr:extended;
              procedure spfl(mask:string;var sf:tspfl; var k:integer);
              function arcsin(x:mtreal):mtreal;
              function arccos(x:mtreal):mtreal;
              function tg(  x:mtreal):mtreal ;
              procedure mappend(var rs:text);







implementation

function mjdr(day,month,year:word;hour:extended):extended;
var a:extended; b:integer;
begin
      a:=10000.0*year+100.0*month+day;
      if month<=2 then begin month:=month+12 ; year:=year-1; end;
      if a<=15821004.1
               then
               begin
               b:=-2+trunc((year+4716)/4)-1179
               end
               else b:=trunc(year/400)-trunc(year/100)+trunc(year/4);
      a:=365.0*year-679004.0;
      mjdr:=a+b+trunc(30.6001*(month+1))+day+hour/24.0;
end;

function tekmjdr:extended;
var Year, Month, Day, DayOfWeek: Word;
 Hour, Minute, Second, Sec100: Word;
 exhour:extended;
begin
     GetDate(Year, Month, Day, DayOfWeek);
     GetTime(Hour, Minute, Second, Sec100);
     exhour:=hour+Minute/60.0+Second/3600.0+Sec100/360000.0;
     tekmjdr:=mjdr(day,month,year,exhour);
end;
procedure spfl(mask:string;var sf:tspfl; var k:integer);
var de:integer;
    hu:searchrec;
begin
      k:=0;
      findfirst(mask,$3F,hu);
      de:=doserror;
      if de=0 then
      begin
           k:=1;
           sf^[1]:=hu.name;
           {ahu[1]:=hu;}
      end;
      while (de=0) and (k<5000) do
      begin
           findnext(hu);
           de:=doserror;
           if de=0 then
           begin
                k:=k+1;
                sf^[k]:=hu.name;
                {ahu[k]:=hu;}
           end;
     end;
end;
function arcsin(x:mtreal):mtreal;
begin
     if (1-x*x)<=0 then arcsin:=pi/2
                      else arcsin:=arctan(x/sqrt(1-x*X));
end;



function arccos(x:mtreal):mtreal;
begin
     {if x=0 then arccos:=pi/2
            else arccos:= arctan(sqrt(1-x*x)/x);}
     if (abs(x-1)<1e-18) or (abs(x+1)<1e-18) then
     begin
          arccos:=0;
          if abs(x+1)<1e-18 then arccos:=pi;
     end
            else arccos:=(pi/2)-arctan(x/sqrt((1.0-x)*(1.0+x)));

end;
function tg(  x:mtreal):mtreal ;
begin
     tg:=sin(x)/cos(x);
end;

procedure mappend(var rs:text);
begin
{$I-}
     append(rs);
     if ioresult<>0 then rewrite(rs);
{$I+}
end;




begin
       oprza:={true}false;

end.