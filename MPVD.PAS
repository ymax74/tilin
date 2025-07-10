unit mpvd ;
interface
         uses sbop;
         procedure pvdv(x,y,z,xt,yt,zt,kp:extended;var a,e,i,omb,omm,w,tp:extended);
         function skal(a,b:vector):extended;
         procedure vecp(a,b:vector;var c:vector);
         procedure sort(spfl:dst2;i:integer);
         procedure sort3(spfl:ttrs;i:integer);
         procedure sort2(spfl:dst2;ai:idst2;i:integer);
         procedure sort6(spfl,aex,aex2:dst2;aex3:dst3;i:integer);
         procedure sort7(spfl,aex,aex2,aex3,aex4:ats;i:integer);

         procedure clan(m1,m2:extended;teta:dst3;x,y,ro:dst3;nnt:integer
         ;ain:idst2;aex,aex2:dst2;gbop,gbop2,gbop3:idst3;aex3:dst3);
         procedure filtrpol(pst:string;nnt:integer;gbop:idst3;x,y,ro:dst3;ain:idst2;var mxk2:integer);
         procedure telclan(y:dst3;nnt:integer;ain:idst2;aex,aex2:dst2;gbop,gbop2,gbop3:idst3);
         procedure telclan3(y:dst3;nnt:integer;ain:idst2;aex,aex2:dst2;gbop,gbop2,gbop3:idst3;t5:extended);
         var tvy:integer;
implementation



procedure sort2(spfl:dst2;ai:idst2;i:integer);
var ii,jj,xi:integer;
x:extended;
me:byte;
hy:boolean;

begin
     for ii:=2 to i do
     begin
          for jj:=i downto ii do
          begin
               hy:=(spfl^[jj-1]>spfl^[jj]);

               if hy then
               begin
                     x:=spfl^[jj-1];
                     spfl^[jj-1]:=spfl^[jj];
                     spfl^[jj]:=x;

                     xi:=ai^[jj-1];
                     ai^[jj-1]:=ai^[jj];
                     ai^[jj]:=xi;
               end;
          end;
      end;
 end;
procedure sort(spfl:dst2;i:integer);
var ii,jj:integer;
x:extended;
me:byte;
hy:boolean;
begin
     for ii:=2 to i do
     begin
          for jj:=i downto ii do
          begin
               hy:=(spfl^[jj-1]>spfl^[jj]);

               if hy then
               begin
                     x:=spfl^[jj-1];
                     spfl^[jj-1]:=spfl^[jj];
                     spfl^[jj]:=x;
               end;
          end;
      end;
 end;
procedure sort3;
var ii,jj:integer;
x:extended;
me:byte;
hy:boolean;
begin
     for ii:=2 to i do
     begin
          for jj:=i downto ii do
          begin
               hy:=(spfl[jj-1]>spfl[jj]);

               if hy then
               begin
                     x:=spfl[jj-1];
                     spfl[jj-1]:=spfl[jj];
                     spfl[jj]:=x;
               end;
          end;
      end;
 end;

procedure sort6;
var ii,jj:integer;
x:extended;
me:byte;
hy:boolean;
begin
     for ii:=2 to i do
     begin
          for jj:=i downto ii do
          begin
               hy:=(spfl^[jj-1]>spfl^[jj]);

               if hy then
               begin
                     x:=spfl^[jj-1];
                     spfl^[jj-1]:=spfl^[jj];
                     spfl^[jj]:=x;

                     x:=aex^[jj-1];
                     aex^[jj-1]:=aex^[jj];
                     aex^[jj]:=x;

                     x:=aex2^[jj-1];
                     aex2^[jj-1]:=aex2^[jj];
                     aex2^[jj]:=x;

                     x:=aex3^[jj-1];
                     aex3^[jj-1]:=aex3^[jj];
                     aex3^[jj]:=x;



               end;
          end;
      end;
 end;

procedure sort7;
var ii,jj:integer;
x:extended;
me:byte;
hy:boolean;
begin
     for ii:=2 to i do
     begin
          for jj:=i downto ii do
          begin
               hy:=(spfl^[jj-1]>spfl^[jj]);

               if hy then
               begin
                     x:=spfl^[jj-1];
                     spfl^[jj-1]:=spfl^[jj];
                     spfl^[jj]:=x;

                     x:=aex^[jj-1];
                     aex^[jj-1]:=aex^[jj];
                     aex^[jj]:=x;

                     x:=aex2^[jj-1];
                     aex2^[jj-1]:=aex2^[jj];
                     aex2^[jj]:=x;

                     x:=aex3^[jj-1];
                     aex3^[jj-1]:=aex3^[jj];
                     aex3^[jj]:=x;

                     x:=aex4^[jj-1];
                     aex4^[jj-1]:=aex4^[jj];
                     aex4^[jj]:=x;



               end;
          end;
      end;
 end;


function skal ;
begin
     skal:=a[1]*b[1]+a[2]*b[2]+a[3]*b[3];
end;

procedure vecp;
begin
      c[1]:=  a[2]*b[3]-b[2]*a[3];
      c[2]:= (a[3]*b[1]-a[1]*b[3]);
      c[3]:= (a[1]*b[2]-b[1]*a[2]);
end;
function por(x,p:extended):extended;
begin
    {writeln('x=',x:10:5,' p=',p:10:5);}

     por:=exp(p*ln(x));

end;



procedure pvdv;
var r,v,rt,ecosw,esinw,cosw,sinw,gamz,cosi,sini,cosomm,sinomm,sqp,dfte,dxw,dyw:extended;
    rb,vb,qq,zz,ombv,nn,pp,hj:vector;

begin
     tvy:=tvy+1;
     (*writeln(tvy{,' ',x,' ',y,' ',z,' ',xt,' ',yt,' ',zt,' ',kp});*)
     if tvy=16 then
     begin
          tvy:=16;
     end;
     r:=sqrt(x*X+y*y+z*z);
     v:=sqrt(xt*xt+yt*yt+zt*zt);

     rb[1]:=x/r;
     rb[2]:=y/r;
     rb[3]:=z/r;

     vb[1]:=xt/v;
     vb[2]:=yt/v;
     vb[3]:=zt/v;

     sini:=skal(rb,vb);
     gamz:=sqrt(1-sini*sini);
     vecp(rb,vb,qq);
     qq[1]:=qq[1]/gamz;
     qq[2]:=qq[2]/gamz;
     qq[3]:=qq[3]/gamz;
     {writeln('Q1=',qq[1]:10:7);
     writeln('Q2=',qq[2]:10:7);
     writeln('Q3=',qq[3]:10:7);

     writeln('mqq=',sqrt(qq[1]*qq[1]+qq[2]*qq[2]+qq[3]*qq[3]):20:10);}
     {readln;}

     a:=kp*kp/( ((2*kp*kp)/r) - v*v);

     vecp(rb,vb,hj);
     {writeln('mrb=',sqrt(rb[1]*rb[1]+rb[2]*rb[2]+rb[3]*rb[3]):20:10);
     writeln('mvb=',sqrt(vb[1]*vb[1]+vb[2]*vb[2]+vb[3]*vb[3]):20:10);}

     if (qq[1])<>0 then
     begin
          sqp:=hj[1]/qq[1];
     end else
     begin
           if (qq[2])<>0 then
           begin
                sqp:=hj[2]/qq[2];
           end else
           begin
                if (qq[3])<>0 then
                begin
                     sqp:=hj[3]/qq[3];
                end else
                begin
                     sqp:=0;
                end;

           end;

     end;
     sqp:=sqp*(r*v)/(kp);
     rt:=v*skal(rb,vb);
     ecosw:=-1+sqp*sqp/r;
     esinw:=sqp*rt/(kp);
     e:=sqrt(ecosw*ecosw+esinw*esinw);
     cosw:=ecosw/e;
     sinw:=esinw/e;
     w:=arccos(cosw);
     if sinw<0 then w:=2*pi-w;

     zz[1]:=0;
     zz[2]:=0;
     zz[3]:=1;
     cosi:=skal(zz,qq);
     (*if y=2.90663593762743E+0001 then
     begin

          {cosi:=1;}
     end;
     writeln('cosi=',cosi:10:5,' y=',y);*)
     if tvy=11 then
     begin
          tvy:=11;
     end;

     i:=arccos((cosi));

     sini:=abs(sin(i));

     vecp(zz,qq,ombv);
     if sini<>0 then
     begin
          ombv[1]:=ombv[1]/sini;
          ombv[2]:=ombv[2]/sini;
          ombv[3]:=ombv[3]/sini;
    end else
    begin
          ombv[1]:=0;
          ombv[2]:=0;
          ombv[3]:=0;

    end;

     dxw:=ombv[1];
     dyw:=ombv[2];
     dfte:=sqrt(sqr(dyw)+sqr((dxw)));
     if dfte<>0 then
     begin
          omb:=arctan(dxw/dyw);
          if dyw<0 then omb:=pi+omb;
          if omb<0 then omb:=omb+2*pi;
     end else
     begin
          omb:=0;
     end;


     vecp(rb,qq,nn);
     pp[1]:=rb[1]*cosw+nn[1]*sinw;
     pp[2]:=rb[2]*cosw+nn[2]*sinw;
     pp[3]:=rb[3]*cosw+nn[3]*sinw;
     {writeln('pp1=',pp[1]:10:5);
     writeln('pp2=',pp[2]:10:5);
     writeln('pp3=',pp[3]:10:5);
     readln;}
     cosomm:=skal(ombv,pp);
     if qq[1]<>0 then
     begin
          sinomm:=(ombv[2]*pp[3]-ombv[3]*pp[2])/qq[1];
     end else
     begin
          sinomm:=0;
     end;
     {writeln('siomm=',sinomm:10:5,180-arcsin(sinomm)*180/pi:10:5);}

     omm:=arcsin(sinomm);
     if cosomm<0 then omm:=pi-omm;
     w:=2*arctan(sqrt(abs((1-e)/(1+e)))*tg(w/2));

     tp:=-(w-e*sin(w))/(kp*por(abs(a),-3/2));
     a:=abs(a);

end;

procedure telclan(y:dst3;nnt:integer;ain:idst2;aex,aex2:dst2;gbop,gbop2,gbop3:idst3);
var j,p2,ncl,mi,ti,i,po1,po2:integer;
grn,t5,raz,mxg:extended;
afl:array[1..300] of integer;

begin
     t5:=20;
     for i:=1 to 300 do afl[i]:=0;
     for i:=1 to nnt do
     begin
          ain^[i]:=i;
     end;

     for i:=1 to nnt do
     begin
          aex^[i]:=y^[i];
     end;
     sort2(aex,ain,nnt);

     for i:=2 to nnt do
     begin
          aex2^[i-1]:=aex^[i]-aex^[i-1];
     end;
     sort(aex2,nnt-1);
     p2:=trunc((nnt-1)*(1-intver));
     p2:=trunc((nnt-1)-(nnt-1)*0.08);
     {writeln('p2=',p2,' nnt=',nnt);}
     grn:=t5*aex2^[p2];
     ncl:=1;
     for i:=1 to nnt do
     begin
          gbop2^[i]:=0;
     end;
     gbop2^[1]:=1;
     afl[1]:=1;
     for i:=2 to nnt do
     begin
          raz:=aex^[i]-aex^[i-1];
          if raz>grn then
          begin
               ncl:=ncl+1;
          end;
          gbop2^[i]:=ncl;
          if ncl<=300 then
          begin
               afl[ncl]:=afl[ncl]+1;
          end;
     end;
     {writeln('ncl=',ncl);}
     if ncl>300 then ncl:=300;
     mxg:=-10000;
     mi:=-1;

     for i:=1 to ncl do
     begin
          if afl[i]>mxg then
          begin
               mxg:=afl[i];
               mi:=i;
          end;
     end;
     po1:=0;
     po2:=0;
     for i:=1 to nnt do
     begin
          ti:=ain^[i];
          if (gbop2^[i]=mi) and (gbop^[ti]=0) then
          begin
               gbop^[ti]:=0;
               po1:=po1+1;
          end else
          begin
               gbop^[ti]:=1;
               po2:=po2+1;
          end;
     end;
     {writeln(po1,' ',po2);}
end;
procedure telclan3;
var j,p2,ncl,mi,ti,i,po1,po2:integer;
grn,raz,mxg:extended;
afl:array[1..300] of integer;

begin
     for i:=1 to 300 do afl[i]:=0;
     for i:=1 to nnt do
     begin
          ain^[i]:=i;
     end;

     for i:=1 to nnt do
     begin
          aex^[i]:=y^[i];
     end;
     sort2(aex,ain,nnt);

     for i:=2 to nnt do
     begin
          aex2^[i-1]:=aex^[i]-aex^[i-1];
     end;
     sort(aex2,nnt-1);
     p2:=trunc((nnt-1)*(1-intver));
     p2:=trunc((nnt-1)-(nnt-1)*0.08);
     {writeln('p2=',p2,' nnt=',nnt);}
     grn:=t5*aex2^[p2];
     ncl:=1;
     for i:=1 to nnt do
     begin
          gbop2^[i]:=0;
     end;
     gbop2^[1]:=1;
     afl[1]:=1;
     for i:=2 to nnt do
     begin
          raz:=aex^[i]-aex^[i-1];
          if raz>grn then
          begin
               ncl:=ncl+1;
          end;
          gbop2^[i]:=ncl;
          if ncl<=300 then
          begin
               afl[ncl]:=afl[ncl]+1;
          end;
     end;
     {writeln('ncl=',ncl);}
     if ncl>300 then ncl:=300;
     mxg:=-10000;
     mi:=-1;

     for i:=1 to ncl do
     begin
          if afl[i]>mxg then
          begin
               mxg:=afl[i];
               mi:=i;
          end;
     end;
     po1:=0;
     po2:=0;
     for i:=1 to nnt do
     begin
          ti:=ain^[i];
          if (gbop2^[i]=mi) and (gbop^[ti]=0) then
          begin
               gbop^[ti]:=0;
               po1:=po1+1;
          end else
          begin
               gbop^[ti]:=1;
               po2:=po2+1;
          end;
     end;
     {writeln(po1,' ',po2);}
end;
procedure tel2clan(y:dst3;nnt:integer;ain:idst2;aex,aex2:dst2;gbop,gbop2,gbop3:idst3;aex3:dst3);
var p1,p2,rp,km,n2,i,j:integer;
v1:extended;
begin
     for i:=1 to nnt do
     begin
            p1:=i-5;
            if p1<1 then p1:=1;
            p2:=p1+10;
            if p2>nnt then
            begin
                 p2:=nnt;
                 p1:=p1-10;
                 if p1<1 then p1:=1;
            end;
            n2:=trunc((p2-p1)/2)+1;
            km:=0;
            for j:=p1 to p2 do
            begin
                 km:=km+1;
                 aex^[km]:=y^[j];
            end;
            sort(aex,km);
            v1:=aex^[n2];
            aex3^[i]:=y^[i]-v1;
     end;
     telclan(aex3,nnt,ain,aex,aex2,gbop,gbop2,gbop3);


end;
procedure clan2(x,y,ro:dst3;nnt:integer;ain:idst2;aex,aex2:dst2;gbop,gbop2,gbop3:idst3;aex3:dst3);
begin
     {exit;}
     telclan(y,nnt,ain,aex,aex2,gbop,gbop2,gbop3);
     telclan(x,nnt,ain,aex,aex2,gbop,gbop2,gbop3);
     telclan(ro,nnt,ain,aex,aex2,gbop,gbop2,gbop3);
     {telclan(et);}
     tel2clan(x,nnt,ain,aex,aex2,gbop,gbop2,gbop3,aex3);

     {tel2clan(x,nnt,ain,aex,aex2,gbop,gbop2,gbop3,aex3);
     tel2clan(x,nnt,ain,aex,aex2,gbop,gbop2,gbop3,aex3);}
     {tel2clan(y,nnt,ain,aex,aex2,gbop,gbop2,gbop3,aex3);
     tel2clan(ro,nnt,ain,aex,aex2,gbop,gbop2,gbop3,aex3);}
end;

procedure clan(m1,m2:extended;teta:dst3;x,y,ro:dst3;nnt:integer;ain:idst2;aex,aex2:dst2;gbop,gbop2,gbop3:idst3;aex3:dst3);
var pdt,pdx,pdy:extended;
j,ib:integer;
begin
     for ib:=1 to nnt do
     begin
          gbop^[ib]:=0;
     end;
     clan2(x,y,ro,nnt,ain,aex,aex2,gbop,gbop2,gbop3,aex3);
     for ib:=1 to nnt do
     begin
          gbop3^[ib]:=gbop^[ib];
     end;

     {writeln(rsv);
     for ib:=1 to nnt do
     begin
          gbop3^[ib]:=gbop^[ib];
          writeln(rsv,gbop^[ib],' ',et^[ib]+t0:8:3,' ',teta^[ib]:10:5,' ',ro^[ib]:10:5,' ',x^[ib]:10:5,' ',y^[ib]:10:5);
     end;
     writeln(rsv);
     flush(rsv);}
     if (abs(m2-m1)<1)  and false then
     begin

          for ib:=1 to nnt do
          begin
               if gbop3^[ib]=1 then
               begin
                    pdt:=teta^[ib];
                    pdx:=x^[ib];
                    pdy:=y^[ib];
                    teta^[ib]:=teta^[ib]-180;
                    if teta^[ib]<0 then teta^[ib]:=teta^[ib]+360;
                    x^[ib]:=ro^[ib]*sin(teta^[ib]*pi/180);
                    y^[ib]:=ro^[ib]*cos(teta^[ib]*pi/180);
                    for j:=1 to nnt do
                    begin
                         gbop^[j]:=0;
                    end;
                    clan2(x,y,ro,nnt,ain,aex,aex2,gbop,gbop2,gbop3,aex3);;

                    if gbop^[ib]=1 then
                    begin
                         teta^[ib]:=pdt;
                         x^[ib]:=pdx;
                         y^[ib]:=pdy;
                    end;

               end;
          end;
     end;
     {writeln(rsv);
     for ib:=1 to nnt do
     begin
          gbop3^[ib]:=gbop^[ib];
          writeln(rsv,gbop^[ib],' ',et^[ib]+t0:8:3,' ',teta^[ib]:10:5,' ',ro^[ib]:10:5,' ',x^[ib]:10:5,' ',y^[ib]:10:5);
     end;
     writeln(rsv);
     flush(rsv);}


end;
procedure filtrpol(pst:string;nnt:integer;gbop:idst3;x,y,ro:dst3;ain:idst2;var mxk2:integer);
var aa,bb,h,x1,y1,x2,y2,r1,r2,teh,mxteh:extended;
mxj,j,n1,n2,i:integer;
mxn,mxn2,kt,irz,mirz,k1,k2:integer;


begin
     aa:=0;
     bb:=360;
     h:=1;
     mxj:=trunc((bb-aa)/h);
     mxn:=0;
     for j:=0 to mxj do
     begin
          teh:=aa+h*j;
          x1:=1000*sin(teh*pi/180);
          y1:=1000*cos(teh*pi/180);

          x2:=1000*sin((teh+180)*pi/180);
          y2:=1000*cos((teh+180)*pi/180);
          n1:=0;
          n2:=0;

          for i:=1 to nnt do
          begin
               if gbop^[i]=0 then
               begin
                    r1:=sqrt(sqr(x^[i]-x1)+sqr(y^[i]-y1));
                    r2:=sqrt(sqr(x^[i]-x2)+sqr(y^[i]-y2));
                    if r1>r2 then
                    begin
                         n1:=n1+1;
                    end else
                    begin
                         n2:=n2+1;
                    end;
               end;

          end;

          if mxn<n1 then
          begin
               mxn:=n1;
               mxteh:=teh;
               mxn2:=n2;
          end;
     end;
          teh:=mxteh;
          x1:=1000*sin(teh*pi/180);
          y1:=1000*cos(teh*pi/180);

          x2:=1000*sin((teh+180)*pi/180);
          y2:=1000*cos((teh+180)*pi/180);
          n1:=0;
          n2:=0;
          kt:=1;
          ain^[kt]:=0;

          for i:=1 to nnt do
          begin
               if gbop^[i]=0 then
               begin
                    r1:=sqrt(sqr(x^[i]-x1)+sqr(y^[i]-y1));
                    r2:=sqrt(sqr(x^[i]-x2)+sqr(y^[i]-y2));
                    if r1>r2 then
                    begin
                         {n1:=n1+1;}
                    end else
                    begin
                         gbop^[i]:=9;
                         kt:=kt+1;
                         ain^[kt]:=i;
                    end;
               end;

          end;
          kt:=kt+1;
          ain^[kt]:=nnt+1;
          if kt>2 then
          begin
              mirz:=0;
              for i:=1 to kt-1 do
              begin
                   irz:=ain^[i+1]-ain^[i];
                   irz:=0;
                   for j:=ain^[i]+1 to ain^[i+1]-1 do
                   begin
                        if (gbop^[j]=0) then irz:=irz+1;
                   end;

                   if mirz<irz then
                   begin
                        mirz:=irz;
                        k1:=ain^[i];
                        k2:=ain^[i+1];
                   end;
              end;
              if mirz>0 then
              begin

              for i:=1 to nnt do
              begin
                    if (gbop^[i]=0) and not( (k1<=i) and (i<=k2)) then
                    begin
                         gbop^[i]:=8;
                    end;
              end;
              end;
          end;
          mxk2:=mxn2;
          {writeln(rsfp,pst,' ',mxn:5,' ',mxn2:5,' ',mxteh:4:0);
          flush(rsfp);}
          {telclan(et);}


     i:=0;


end;


begin
     tvy:=0;
end.

