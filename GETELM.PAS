unit getelm;
interface
          procedure getel(swds,comp10:string;var abp,ev,th,prd,ai,omm,omb,ro2014,teta2014:extended);
implementation
uses dos;
const otm='..\wds';
      t600=600;

type pkt=array[1..t600] of string[12];
     tspfl=^pkt;

var hj:integer;
sf4,awds:tspfl;
k4:integer;
j2000:extended;
uwds:boolean;


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

procedure  crtspwds;
var i,j,mi:integer;
dt:text;
pst30:string[30];
pst10:string[10];

begin
     spfl(otm+'\w*.txt',sf4,k4);
     for i:=1 to k4 do
     begin
          assign(dt,otm+'\'+sf4^[i]);
          reset(dt);
          readln(dt,pst30,pst10);
          awds^[i]:=pst10;
          close(dt);
     end;
     i:=0;
end;
procedure erval(fed:string; ercode:integer);
begin
     if ercode>0 then
     begin
          writeln(' stroka=',fed,' code=',ercode);
          {readln;}
     end;
end;
function otnc(rf:string):string;
var pst:string;
   k,i,l:integer;
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
     l:=length(pst);
     if pst[l]='.' then
     begin
          pst:=copy(pst,1,l-1);
          l:=0;
     end;
     otnc:=pst;
end;

procedure getel(swds,comp10:string;var abp,ev,th,prd,ai,omm,omb,ro2014,teta2014:extended);
var mi,i,j,k,er:integer;
dt:text;
pst,pst2:string;
ch:char;
comp5:string[5];

begin
     if not(uwds) then
     begin
          crtspwds;
          uwds:=true;
     end;

{
ORBITAL ELEMENTS:

                  ------P-------  ---a------  --i-----  Node----  -----T-------  -e------  omega---  G  Reference
HLD  60              223.2     y    0.8798 a  128.3     324.7      1903.08    y  0.6479    325.7     3  Hrt2011a
               +/-    12.2          0.0039      2.8       2.2         1.28       0.0123      3.8

                                       2014.0           2015.0           2016.0           2017.0
STF1540AB         5    Hop1960a    146.9  28.632    146.9  28.625    146.9  28.619    146.8  28.612
}
     comp5:=copy(comp10,1,5);
     mi:=-1;
     for i:=1 to k4 do
     begin
          if swds=awds^[i] then mi:=i;
     end;
     abp:=-1;
     if mi>0 then
     begin
          assign(dt,otm+'\'+sf4^[mi]);
          reset(dt);
          while not(eof(dt)) do
          begin
               pst:='';
               ch:=' ';
               while ch<>#10 do
               begin
                    read(dt,ch);
                    pst:=pst+ch;
               end;
               if      (pst[1]<>' ')
                  and  (pst[39]='.')
                  and  (pst[38]<>'.')
                  and  (pst[45]='.')
                  and  (pst[56]='.')
                  and  (pst[62]='.')
                  then
               begin
                    pst2:=copy(pst,8,5);
                    if comp5=pst2 then
                    begin
                    pst2:=copy(pst,36,41-36+1);
                    pst2:=otnc(pst2);
                    val(pst2,teta2014,er);
                    erval(pst2,er);
                    pst2:=copy(pst,42,48-42+1);
                    pst2:=otnc(pst2);
                    val(pst2,ro2014,er);
                    erval(pst2,er);
                    end;
               end;

               if      (pst[1]<>' ')
                  and  (pst[25]='.')
                  and  (pst[38]='.')
                  and  (pst[50]='.')
                  and  (pst[60]='.')
                  and  (pst[72]='.')
                  and  (pst[83]='.')
                  and  (pst[95]='.')
                  then
               begin
                    pst2:=copy(pst,8,5);
                    if comp5=pst2 then
                    begin
                    pst2:=copy(pst,19,31-19+1);
                    pst2:=otnc(pst2);
                    val(pst2,prd,er);
                    erval(pst2,er);
                    if pst[32]='d' then
                    begin
                         if prd=8125.4 then prd:=prd+10000;
                         prd:=prd/365.25;
                    end;
                    if pst[32]='c' then
                    begin
                         prd:=prd*100;
                    end;

                    pst2:=copy(pst,35,43-35+1);
                    pst2:=otnc(pst2);
                    val(pst2,abp,er);
                    erval(pst2,er);

                    pst2:=copy(pst,47,54-47+1);
                    pst2:=otnc(pst2);
                    val(pst2,ai,er);
                    erval(pst2,er);

                    pst2:=copy(pst,57,64-57+1);
                    pst2:=otnc(pst2);
                    val(pst2,omb,er);
                    erval(pst2,er);

                    pst2:=copy(pst,67,78-67+1);
                    pst2:=otnc(pst2);
                    val(pst2,th,er);
                    erval(pst2,er);
                    if pst[79]='m' then
                    begin
                         th:=th-j2000;
                         th:=2000+th/365.25;
                         er:=0;
                    end;

                    pst2:=copy(pst,82,89-82+1);
                    pst2:=otnc(pst2);
                    val(pst2,ev,er);
                    erval(pst2,er);

                    pst2:=copy(pst,92,99-92+1);
                    pst2:=otnc(pst2);
                    val(pst2,omm,er);
                    erval(pst2,er);
                    end;
                    i:=23;
               end;

          end;
          close(dt);

     end;

end;
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


begin
     new(sf4);new(awds);
     j2000:=mjdr(1,1,2000,12);
     {crtspwds;}
     uwds:=false;
end.
