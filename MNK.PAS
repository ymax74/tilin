{$A+,B-,D+,E-,F+,G+,I+,L+,N+,O-,P-,Q+,R+,S+,T-,V+,X+,Y+}
{$M 16384,0,655360}
unit mnk;
interface
          const
               {nk2=16000;}
                nk2=1500;
                {nk2=500;}

               s=7;
               {nk=1240;}
               {nk=2000;}
               nk=2000;
               {nk=500;}
          type
               matrls=array[1..s+3,1..s+1+3] of extended;
               stb=array[1..s+3] of extended;
               stb2=array[0..s+3] of extended;
               mslss=array[1..nk,1..s+1] of single;
               msls=^mslss;
               mslss2=array[1..s+1] of extended;
               msls20=array[1..nk2] of ^mslss2;
               msls2=^msls20;
                 dft=array[1..22,1..375] of word;
                 juy=^dft;
                 trh=array[1..11] of juy;
              trdf= array[0..s+3,0..s+1+3] of extended;

               {mslss2=array[1..s+1] of single;
               mslss3=array[1..nk] of ^mslss2;
               msls2=^mslss3;}
                tr={procedure}function(i,j:integer):extended;
          var km0:integer;
          procedure  slvls(var b:matrls;s:integer;var x:stb);
          procedure  slvls2(var b:matrls;s:integer;var x:stb);
          procedure  slvls3(var b:matrls;n:integer;var x2:stb);

          procedure  ptlst(a: msls;var x:stb;var c:matrls;km,sn:integer);
          procedure  ptlst4(a: msls2;var x:stb;var c:matrls;km,sn:integer);

          procedure  ptlst2(a: msls2;var x:stb;var c:matrls;km,sn:integer);
          procedure  ptlst3(a: msls;var x:stb;var c:matrls;km,sn:integer;var c3:matrls);


          Function Ranorm(s:extended):extended;
          function det(var b:matrls;s:integer):extended;
          procedure sgmnk(a:msls;c:matrls;x:stb;var sigm:stb; sn,km:integer);
          procedure sgmnk2(a:msls2;c:matrls;x:stb;var sigm:stb; sn,km:integer);

          {procedure sgmnk2(a:msls2;c:matrls;x:stb;var sigm:stb; sn,km:integer);}


          var alfa,beta,gama,rz,zk,xcet,ycet,
    rk,kz,sdfx,sdfy,sdfz,yp,ik,jk,bb,as,bs,bbs,z,ki,bkl,w:extended;
   i,xc,yc,y2,y1,j,km,k,xmt,ymt:integer;
   xf2,sgm,xft:stb;
   cv:matrls;
   rze:text;
   {a:msls2;}
   il:longint;
   imgt:trh;
   {appc:tr;}


implementation

function app(i,j:integer):extended;
begin
                     (*ik:=i;
                     jk:=j;
                     km:=km+1;
                     rz:=xf[4]*sqr(ik-xce)+xf[5]*sqr(jk-yce);
                     beta:=exp((alfa-1)*ln(rz));
                     gama:=beta*rz;
                     z:=1/(beta*rz+bb);
                     beta:=alfa*beta;
                     kz:=-(z*z);
                     a^[km,1]:=2*as*sqr(i-xce)*kz*beta;
                     a^[km,2]:=2*bs*sqr(j-yce)*kz*beta;
                     a^[km,3]:=-2*xf[4]*(i-xce)*kz*beta;
                     a^[km,4]:=-2*xf[5]*(j-yce)*kz*beta;
                     a^[km,5]:=2*bbs*1*kz;
                     a^[km,6]:=1{gama*ln(rz)};
                     z:=z+bkl;
                     a^[km,6]:=z-img[(i-1) div 22 +1]^[(i-1) mod 22 +1,j];*)

                     ik:=(i-1) div 242+1;
                     jk:=(i-1) mod 242 +1;
                     rz:=xft[4]*sqr(ik-xcet)+xft[5]*sqr(jk-ycet);
                     beta:=exp((alfa-1)*ln(rz));
                     gama:=beta*rz;
                     z:=1/(beta*rz+bb);
                     beta:=alfa*beta;
                     kz:=-(z*z);
                     if j=1 then app:=2*as*sqr(i-xceT)*kz*beta;
                     if j=2 then app:=2*bs*sqr(j-yceT)*kz*beta;
                     if j=3 then app:=-2*xfT[4]*(i-xceT)*kz*beta;
                     if j=4 then app:=-2*xfT[5]*(j-yceT)*kz*beta;;
                     if j=5 then app:=2*bbs*1*kz;;
                     if j=6 then app:=z-imgt[(i-1) div 22 +1]^[(i-1) mod 22 +1,j];;



end;


function det(var b:matrls;s:integer):extended;
var
   i,k,m,ii,kk:integer;
   a:matrls;
   mul:extended;
begin

      a:=b;
      for m:=1 to s-1 do
      begin
            for i:=m+1 to s do
            begin
                 for k:=m+1 to s+1 do
                 begin
                       a[i,k]:=a[i,k]-a[i,m]*a[m,k]/a[m,m];
                 end;
            end;
      end;
      mul:=1;
      for m:=1 to s do mul:=mul*a[m,m];
      det:=mul;
end;


procedure sgmnk(a:msls;c:matrls;x:stb;var sigm:stb; sn,km:integer);

var sum,sum2,d:extended;
    i,j,k,mn:integer;
    cd:matrls;

begin
     sum:=0;
     for i:=1 to km do
     begin
          sum2:=0;
          for j:=1 to sn do
          begin
               sum2:=sum2+a^[i,j]*x[j];
          end;
          sum:=sum+sqr(a^[i,sn+1]-sum2);
     end;
     {d:=det(c,sn);}
     for j:=1 to sn do
     begin
          {for k:=1 to sn-1 do
          for mn:=1 to sn-1 do
          begin
               if (mn<j) and  (k<j) then cd[k,mn]:=c[k,mn];
               if (mn<j) and  (k>=j) then cd[k,mn]:=c[k+1,mn];
               if (k<j)  and  (mn>=j) then cd[k,mn]:=c[k,mn+1];
               if (k>=j)  and  (mn>=j) then cd[k,mn]:=c[k+1,mn+1];
          end;}
          sigm[j]:=sqrt(abs(sum/((km-sn)*c[j,j])));
          {writeln(sum/(km-sn));
          }
     end;
     sigm[j+1]:=sqrt(sum/((km-sn)));
end;

procedure sgmnk2(a:msls2;c:matrls;x:stb;var sigm:stb; sn,km:integer);

var sum,sum2,d:extended;
    i,j,k,mn:integer;
    cd:matrls;

begin
     sum:=0;
     for i:=1 to km do
     begin
          sum2:=0;
          for j:=1 to sn do
          begin
               sum2:=sum2+a^[i]^[j]*x[j];
          end;
          sum:=sum+sqr(a^[i]^[sn+1]-sum2);
     end;
     {d:=det(c,sn);}
     for j:=1 to sn do
     begin
          {for k:=1 to sn-1 do
          for mn:=1 to sn-1 do
          begin
               if (mn<j) and  (k<j) then cd[k,mn]:=c[k,mn];
               if (mn<j) and  (k>=j) then cd[k,mn]:=c[k+1,mn];
               if (k<j)  and  (mn>=j) then cd[k,mn]:=c[k,mn+1];
               if (k>=j)  and  (mn>=j) then cd[k,mn]:=c[k+1,mn+1];
          end;}
          sigm[j]:=sqrt(sum/((km-sn)*c[j,j]));
          {writeln(sum/(km-sn));
          }
     end;
     sigm[j+1]:=sqrt(sum/((km-sn)));
end;





procedure  slvls(var b:matrls;s:integer;var x:stb);
var
   i,k,m,ii,kk:integer;
   a:matrls;
begin

      a:=b;
      for m:=1 to s-1 do
      begin
            for i:=m+1 to s do
            begin
                 for k:=m+1 to s+1 do
                 begin
                       a[i,k]:=a[i,k]-a[i,m]*a[m,k]/a[m,m];
                 end;
            end;
      end;
      for m:=s downto 1 do
      begin
           x[m]:=a[m,s+1];
           for k:=m+1 to s do x[m]:=x[m]-x[k]*a[m,k];
           x[m]:=x[m]/a[m,m];

      end;
end;

procedure  slvls2(var b:matrls;s:integer;var x:stb);
var
   i,k,m,ii,kk:integer;
   a:matrls;
begin

      a:=b;
      for m:=1 to s-1 do
      begin
            for i:=m+1 to s do
            begin
                 for k:=m+1 to s+1 do
                 begin
                       a[i,k]:=a[i,k]-a[i,m]*a[m,k]/a[m,m];
                 end;
            end;
      end;
      for m:=s downto 1 do
      begin
           x[m]:=a[m,s+1];
           for k:=m+1 to s do x[m]:=x[m]-x[k]*a[m,k];
           x[m]:=x[m]/a[m,m];

      end;
      b:=a;
end;



 Function Ranorm(s:extended):extended;
 var a : extended;
     i : byte;
 begin
  a:=0;
  for i:=1 to 12 do a:=a+random;
  Ranorm:=(a-6.0)*s;
 end;



procedure  ptlst(a: msls;var x:stb;var c:matrls;km,sn:integer);
var
   sum:extended   ;
   i,j,k:integer;
begin
     for k:=1 to sn do
     begin
           for j:=1 to sn+1  do
           begin
                sum:=0;
                for i:=1  to  km  do sum:=sum+a^[i,k]*a^[i,j];
                c[k,j]:=sum;
          end;
     end;
     slvls2(c,sn,x);
end;
procedure  ptlst2;
var
   sum:extended;
   i,j,k:integer;
begin
     for k:=1 to sn do
     begin
           for j:=1 to sn+1  do
           begin
                sum:=0;
                for i:=1  to  km  do sum:=sum+a^[i]^[k]*a^[i]^[j];
                c[k,j]:=sum;
          end;
     end;
     slvls2(c,sn,x);
end;

procedure  ptlst3;
var
   sum:extended   ;
   i,j,k:integer;
begin
     for k:=1 to sn do
     begin
           for j:=1 to sn+1  do
           begin
                sum:=0;
                for i:=1  to  km  do sum:=sum+a^[i,k]*a^[i,j];
                c[k,j]:=sum;
          end;
     end;
     c3:=c;
     slvls2(c,sn,x);
end;
procedure  ptlst4(a: msls2;var x:stb;var c:matrls;km,sn:integer);
var
   sum:extended   ;
   i,j,k:integer;
begin
     for k:=1 to sn do
     begin
           for j:=1 to sn+1  do
           begin
                sum:=0;
                for i:=1  to  km  do sum:=sum+a^[i]^[k]*a^[i]^[j];
                c[k,j]:=sum;
          end;
     end;
     slvls3(c,sn,x);
end;


procedure slvls3(var b:matrls;n:integer;var x2:stb);
{ Variables declaration }
var
  row_with_max_element, i, j, k: integer;
  swp, tmp, max: real;
  matrix:trdf;
  x: stb2;
  ir,jr,n1:integer;
begin
     n1:=n+1;
     for ir:=1 to n do
     begin
          for jr:=1 to n+1 do
          begin
               matrix[ir-1,jr-1]:=b[ir,jr];
          end;
     end;

  {// float max,swp,tmp,x[N1];
  // register short int row_with_max_element,i,j,k;}


  { Printing given matrix }
  {ShowMatrix;}
  { Main loop }
  for i := 0 to pred(N) do begin
    { Searching for max element in the current column (i) }

    max := abs(matrix[i][i]);
    row_with_max_element := i;
    for j := i+1 to pred(N) do
      if abs(matrix[j][i]) > max then begin
        row_with_max_element := j;
        max := abs(matrix[j][i]);
      end;

    { Swapping 2 lines of matrix - row_with_max_element & i }
    if row_with_max_element <> i then
      for j := 0 to pred(N1) do begin
        swp := matrix[row_with_max_element][j];
        matrix[row_with_max_element][j] := matrix[i][j];
        matrix[i][j] := swp;
      end;

    {// Excluding variable x[i] from equations}
    tmp := matrix[i][i];
    for j := N downto i do
      matrix[i][j] := matrix[i][j] / tmp;

    for j := i + 1 to pred(N) do begin
      tmp := matrix[j][i];
      for k := N downto i do
        matrix[j][k] := matrix[j][k] - tmp*matrix[i][k];
    end;
  end;

  { Calculating vector x }
  x[N-1] := matrix[N-1][N];
  for i := N-2 downto 0 do begin
    x[i] := matrix[i][N];
    for j := i + 1 to pred(N) do
      x[i] := x[i] - matrix[i][j]*x[j];
  end;

  { Printing solution }
  {writeln;
  writeln('Solution:');
  for i := 0 to pred(N) do
    writeln('x', (i+1):2, ' = ', x[i]:10:5);}
    for ir:=1 to n do x2[ir]:=x[ir-1];
end;



begin
     km0:=0;
end.
