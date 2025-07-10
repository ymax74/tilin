program preces;
var dt,rs:text;
wds10,pst12,srh:string[10];
rh,rm,rss,dg,dm,ds,ra,dc:extended;
er:integer;
ep,te,ro:extended;
function cyf(rf:char):boolean;
var tc:boolean;
begin
     tc:= (rf='1') or (rf='2') or (rf='3') or (rf='4') or (rf='5') or
     (rf='6') or (rf='7') or (rf='8') or (rf='9') or  (rf='0');
     cyf:=tc;
end;
procedure erval(fed:string; ercode:integer);
begin
     if ercode>0 then
     begin
          writeln(' stroka=',fed,' code=',ercode);
          readln;
     end;

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
     corinp('source.txt');
     assign(dt,'source.txt');
     reset(dt);
     assign(rs,'input.txt');
     rewrite(rs);
     readln(dt,wds10);
     if (length(wds10)<>10) or
     not((wds10[6]='+') or (wds10[6]='-')) or
     not cyf(wds10[1]) or
     not cyf(wds10[2]) or
     not cyf(wds10[3]) or
     not cyf(wds10[4]) or
     not cyf(wds10[5]) or
     not cyf(wds10[7]) or
     not cyf(wds10[8]) or
     not cyf(wds10[9]) or
     not cyf(wds10[10])
     then
     begin
           writeln('Invalid WDS designation');
           readln;
           exit;
     end;
     pst12:=wds10;
     srh:=copy(pst12,1,2);
     val(srh,rh,er);
     erval(srh,er);

     srh:=copy(pst12,3,3);
     val(srh,rm,er);
     erval(srh,er);
     rm:=rm/10.0;

     srh:=copy(pst12,7,2);

     val(srh,dg,er);
     erval(srh,er);

     srh:=copy(pst12,9,2);
     val(srh,ds,er);
     erval(srh,er);

     if pst12[6]='-' then
     begin
          dg:=-dg;
          ds:=-ds;
     end;

     ra:=(rh+rm/60)*15*pi/180;
     dc:=(dg+ds/60)*pi/180;

     ra:=sin(ra);
     dc:=1/cos(dc);



     while not(eof(dt)) do
     begin
          readln(dt,ep,te,ro);
          te:=te+0.00556*ra*(2000-ep)*dc;
          writeln(rs,ep:10:5,' ',te:9:5,' ',ro:9:5);

     end;
     close(dt);
     close(rs);
end.