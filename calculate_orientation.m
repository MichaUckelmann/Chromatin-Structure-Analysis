%variable 'table' needs to be loaded into the workspace and needs to follow
%a Dynamo crop table convention
v = [1; 0; 0]; %defines orientation vector pointing in the direction of the nucleosome face of the template used for template matching
coord = table(:,24:26); %grabs coordinates from dynamo table
[eIdx,eD] = knnsearch(coord,coord,'K',151,'Distance','euclidean'); %finds 150 nearest neighbours for each nucleosome
cart_eul(:,1:3) = table(eIdx(:,1),24:26);
cart_eul(:,4:6) = table(eIdx(:,2),24:26); % select table(eIdx(:, n+1),24:36) for n-th nearest neighbour
cart_eul(:,7:9) = table(eIdx(:,1),7:9);
cart_eul(:,10:12) = table(eIdx(:,2),7:9); %select table(eIdx(:, n+1),7:9) for n-th nearest neighbour
%carthesian coordinates of nearest neighbour in colmuns 1:3, 4:6 of cart_eul
%euler angles of nearest neighburs in columns 7:9, 10:12
cart_eul(:,13:15) = [(cart_eul(:,4)-cart_eul(:,1)),(cart_eul(:,5)-cart_eul(:,2)),(cart_eul(:,6)-cart_eul(:,3))];
%calculates vector from NCP1 pointing towards the nearest neighbour NCP2
%and puts it into column 13:15 in cart_eul
for i = 1:size(cart_eul,1)
    cart_eul(i,16:18) = cart_eul(i,13:15)/norm(cart_eul(i,13:15));
end
%converts vector pointing from NCP1 to NCP 2 to ubit vector and writes it
%into column 16:18 of cart_eul
for i = 1:size(cart_eul,1)
cart_eul(i,19:21) = transpose(dynamo_euler2matrix(cart_eul(i,7:9))*v);
cart_eul(i,22:24)= transpose(dynamo_euler2matrix(cart_eul(i,10:12))*v);
end
%calculates for every NCP1 and NCP2 the rotation matrix using euler angles
%in column 7:9 in cart_eul (NCP1) and in column 10:12 (NCP2). Then applies
%the rotation matrix to a vector v. v should be the unit vector along the axis pointing
%away from the face of the nucleosome. this may change depending on the template used for
%template matching. In this case it is the unit vector along the x-axis v =
%[1; 0; 0]
% cart_eul now contains the unit vectors:
%pointing from NCP1 to NCP2 in columns 16:18
% unit vector pointing into x-direction (NCP-face) of nucleosom 1 in columns
% 19:21
%unit vector pointing in x-direction (NCP-face) of nucleosome 2 in columns
%22:24
cos_alpha = dot(cart_eul(:,19:21),cart_eul(:,22:24),2);
cos_betai = dot(cart_eul(:,16:18),cart_eul(:,19:21),2);
cos_betaj = dot(cart_eul(:,16:18),cart_eul(:,22:24),2);
cart_eul(:,25) = acosd(cos_alpha);
cart_eul(:,26) = acosd(cos_betai);
cart_eul(:,27) = acosd(cos_betaj);
ff_ss = cart_eul(:,25)<45 | cart_eul(:,25)>135;
fs = not(ff_ss);
ff = (cart_eul(:,25)<45 | cart_eul(:,25)>135) & (cart_eul(:,26)<45 | cart_eul(:,26)>135 | cart_eul(:,27)<45 | cart_eul(:,27)>135);
ss = not(ff) & not(fs);
cart_eul(:,28) = fs; %logical for face-to-side NCP interaction, colunm 28
cart_eul(:,29) = ff; % logical for face-to-face NCP interaction, column 29
cart_eul(:,30) = ss; % logical for side-to-side NCP interaction, column 30
second_neighbour = cart_eul;
y2 = [sum(second_neighbour(:,28)),sum(second_neighbour(:,29)),sum(second_neighbour(:,30))];