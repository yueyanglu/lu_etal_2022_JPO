%% plot 


%---------------------------------------- curl & div of uv
 % curl/div of orig flx
 curl = calc_curl_HYCOM(fuE,fvE,scqx,scqy,1);
 div = calc_div_HYCOM(uGS,vGS,scvx,scuy,scpx.*scpy,1);
  div = calc_div_HYCOM(uGS,vGS,ones(JDM,IDM),ones(JDM,IDM),ones(JDM,IDM),1);

 % curl/div of div comp
 curl_phi = calc_curl_HYCOM(uflxD,vflxD,scqx,scqy,1);
 div_phi = calc_div_HYCOM(uflxD,vflxD,scvx,scuy,scpx.*scpy,1);
 
 % curl/div of rot comp
 curl_psi = calc_curl_HYCOM(uflxR,vflxR,scqx,scqy,1);
 div_psi = calc_div_HYCOM(uflxR,vflxR,scvx,scuy,scpx.*scpy,1);

%---------------------------------------- uv onto p-point for plot
% [u_psi_p,v_psi_p] = uv2p(u_psi,v_psi);
% [u_phi_p,v_phi_p] = uv2p(u_phi,v_phi);

%% plot the optimal psi and phie 

cmapstr = 'haline';

%----------------- curl and div
clim = [-5 5];
[psi_plt,phi_plt] = deal(div,phi);
psi_plt(fland) = NaN;
phi_plt(fland) = NaN;

figure
[ha, ~] = tight_subplot(1,2,[.04 .01],[.10 .05],[.05 .02]);
axes(ha(1));
h = pcolor(qlon, qlat, psi_plt); set(h,'EdgeColor', 'none');
set(gca,'xlim',qlon([1 end]),'ylim',qlat([1 end]));
title('\psi','fontsize',14)
cmocean(cmapstr)
caxis(clim)
colorbar
axis equal
axes(ha(2));
h = pcolor(plon, plat, phi_plt); set(h,'EdgeColor', 'none');
set(gca,'xlim',plon([1 end]),'ylim',plat([1 end]));
title('\phi','fontsize',14)
set(h,'EdgeColor', 'none');
cmocean(cmapstr)
caxis(clim)
colorbar
axis equal

set(ha(2),'YTickLabel','');

%% plot the TWO components of flux (u/v_phi, u/v_psi)

cmapstr = 'balance';
% cmapstr = 'deep';
titleStr = {'div of cflx', 'd(ch)/dt'};
titleStr = {'Fu-D', 'Fv-D'};
titleStr = {'Fu-R', 'Fv-R'};

%-------- plot orig flx
clim = [-1e2 1e2];
% clim = [-1e-2 1e-2];
% clim = [0 5e2];
% clim = [-1e-5 1e-5];

% [plt11,plt12] = deal(uflx_al{1},uflx_al{2});  %uch,vch
% [plt11,plt12] = deal(dpm_al{1},dpm_al{2});  %uch,vch
% [plt11,plt12] = deal(dpi_al{1},dpi_al{2});  %uch,vch

% [plt11,plt12] = deal(uflx_zsm1_al{1},vflx_zsm1_al{1});  %uch,vch
% [plt11,plt12] = deal(uGS(:,:,1,end),vGS(:,:,1,end));  %uch,vch
% [plt11,plt12] = deal(lmdu,lmdv);  %uch,vch
% [plt11,plt12] = deal((lmdu - lmdu_iso) ./ lmdu, (lmdv - lmdv_iso) ./ lmdv);  %uch,vch
% [plt11,plt12] = deal(usce,vsce);  %uch,vch
% [plt11,plt12] = deal(div_cflx, dch);  %uch,vch
% [plt11,plt12] = deal(adv_cflx, hdc);  %uch,vch
% [plt11,plt12] = deal(b1, b2);  %uch,vch
% [plt11,plt12] = deal(div, vGS);  %uch,vch
[plt11,plt12] = deal(uflxR,vflxR);
% [plt11,plt12] = deal(uflx,vflx);  %uch,vch
% [plt11,plt12] = deal(uflxD,vflxD);  %uch,vch
% [plt11,plt12] = deal(dpioE, dpioS);  %uch,vch
% [plt11,plt12] = deal((d_hsce) / 0.5/24/3600, (d_hecs)/ 0.5/24/3600);  %uch,vch

% [plt11,plt12] = deal((dpioE - dpiiE) / 0.5/24/3600, (dpioS - dpiiS)/ 0.5/24/3600);  %uch,vch
% [plt11,plt12] = deal(uflxE,vflxE);  %uch,vch
% [plt11,plt12] = deal(dcdx_u(:,:,1),dcdy_v(:,:,1));
win_len = 101;
% plt11 = plt11 - smooth_geom_HYCOM(plt11,scu2,win_len,win_len);
% plt12 = plt12 - smooth_geom_HYCOM(plt12,scv2,win_len,win_len);

% plt11(abs(plt11)<1e-10) = NaN;
% plt12(abs(plt12)<1e-10) = NaN;

figure
[ha, ~] = tight_subplot(1,2,[.04 .01],[.10 .05],[.05 .02]);
axes(ha(1));
plot_field_model(plt11,ulon1d,ulat1d,cmapstr)
title(titleStr{1},'fontsize',14)
caxis(clim)
axes(ha(2));
plot_field_model(plt12,vlon1d,vlat1d,cmapstr)
title(titleStr{2},'fontsize',14)
caxis(clim)
cb = colorbar;
set(cb,'Location','SouthOutside','Position',[0.28 0.22 0.46 0.04])

%%

cmapstr = 'curl';
clim = [-5e-4 5e-4];
% plt11 = dpm_al{1}(:,:,1,1); % div_cflx(:,:,1) - ans;
plt11 = div_cflx; % div_cflx(:,:,1) - ans;

figure
plot_field_model(plt11,plon1d,plat1d,cmapstr)
title('Fu','fontsize',14)
caxis(clim)
colorbar

%%
cmapstr = 'balance';

%-------- plot TWO comp
clim = [-1e1 1e1
        -2e0 2e0];
    
[plt11,plt12,plt21,plt22] = deal(u_psi,v_psi,u_phi,v_phi);
plt11(fland) = NaN;
plt12(fland) = NaN;
plt21(fland) = NaN;
plt22(fland) = NaN;

figure
% 
[ha, ~] = tight_subplot(2,2,[.03 .01],[.03 .05],[.02 .07]);
axes(ha(1));
plot_field_model(plt11,ulon1d,ulat1d,cmapstr)
title('Fu_{\psi}','fontsize',14)
caxis(clim(1,:))
% 
axes(ha(2));
plot_field_model(plt12,vlon1d,vlat1d,cmapstr)
title('Fv_{\psi}','fontsize',14)
caxis(clim(1,:))
cb = colorbar;
set(cb,'Location','EastOutside','Position',[0.94 0.55 0.017 0.33])
% 
axes(ha(3));
plot_field_model(plt21,ulon1d,ulat1d,cmapstr)
title('Fu_{\phi}','fontsize',14)
caxis(clim(2,:))
% 
axes(ha(4));
plot_field_model(plt22,vlon1d,vlat1d,cmapstr)
title('Fv_{\phi}','fontsize',14)
caxis(clim(2,:))
cb = colorbar;
set(cb,'Location','EastOutside','Position',[0.94 0.10 0.017 0.33])

%% plot the curl and div of components of the flux 

cmapstr = 'amp';
titlestr = ["err in u-flx" , "err in v-flx"];

%----------------- curl and div
clim = [0 1e-4];
[curl_plt,div_plt] = deal(err_u,err_v);


figure
[ha, ~] = tight_subplot(1,2,[.04 .01],[.10 .05],[.05 .02]);
axes(ha(1));
h = pcolor(qlon, qlat, curl_plt); set(h,'EdgeColor', 'none');
set(gca,'xlim',qlon([1 end]),'ylim',qlat([1 end]));
title(char(titlestr(1)),'fontsize',14)
cmocean(cmapstr)
caxis(clim)
colorbar
axis equal
axes(ha(2));
h = pcolor(plon, plat, div_plt); set(h,'EdgeColor', 'none');
set(gca,'xlim',plon([1 end]),'ylim',plat([1 end]));
title(char(titlestr(2)),'fontsize',14)
set(h,'EdgeColor', 'none');
cmocean(cmapstr)
caxis(clim)
colorbar
axis equal


%% plot errors in curl & div (2-by-2 plots)

cmapstr = ["amp" , "amp"];
clim = {[0 1e-5],[0 1e-5],[0 1e-2],[0 1e-5]};
% [curl_plt,div_plt,curl_df_plt,div_df_plt] = deal(curl_phi,div_psi,...
%     curl_psi, div_phi);
% [curl_plt,div_plt,curl_df_plt,div_df_plt] = deal(curl_phi,div_psi,...
%     abs((curl_psi - curl) ./ curl), abs((div_phi - div) ./ div));
[curl_plt,div_plt,curl_df_plt,div_df_plt] = deal(curl_phi,div_psi,...
    err_curl_psi, err_div_phi);

figure
% 
subplot(221)
h = pcolor(qlon, qlat, curl_plt); set(h,'EdgeColor', 'none');
set(gca,'xlim',qlon([1 end]),'ylim',qlat([1 end]));
title('curl of F_{\phi}','fontsize',14)
cmocean(char(cmapstr(1)))
caxis(clim{1})
colorbar
axis equal
% 
subplot(222)
h = pcolor(plon, plat, div_plt); set(h,'EdgeColor', 'none');
set(gca,'xlim',plon([1 end]),'ylim',plat([1 end]));
title('div of F_{\psi}','fontsize',14)
cmocean(char(cmapstr(1)))
caxis(clim{2})
colorbar
axis equal
% 
subplot(223)
h = pcolor(qlon, qlat, curl_df_plt); set(h,'EdgeColor', 'none');
set(gca,'xlim',qlon([1 end]),'ylim',qlat([1 end]));
title('Err in curl of F_{\psi}','fontsize',14)
cmocean(char(cmapstr(2)))
caxis(clim{3})
colorbar
axis equal
% 
subplot(224)
h = pcolor(plon, plat, div_df_plt); set(h,'EdgeColor', 'none');
set(gca,'xlim',plon([1 end]),'ylim',plat([1 end]));
title('Err in div of F_{\phi}','fontsize',14)
cmocean(char(cmapstr(2)))
caxis(clim{4})
colorbar
axis equal

%%
jj = 1:2:JDM-1;
ii = 1:2:IDM-1;

figure
subplot(121)
quiver(plon(jj,ii), plat(jj,ii), u_psi_p(jj,ii),v_psi_p(jj,ii),'k');
set(gca,'xlim',plon([1 end]),'ylim',plat([1 end]));
title('vel_{\psi}')
axis equal
subplot(122)
quiver(plon(jj,ii), plat(jj,ii), u_phi_p(jj,ii),v_phi_p(jj,ii),'k');
set(gca,'xlim',plon([1 end]),'ylim',plat([1 end]));
title('vel_{\phi}')
axis equal

figure
quiver(plon(jj,ii), plat(jj,ii), u(jj,ii),v(jj,ii),'k');
set(gca,'xlim',plon([1 end]),'ylim',plat([1 end]));
axis equal

%%


