% ============================================================
% Fig 1: Currents and PCC Voltage
% ============================================================
if update_plots(1)
    figure(1); clf; set(gcf,'WindowStyle','docked');

    subplot(3,1,1)
    plot(t,x1(state_idx.i1d,:),'-', ...
         t,x1(state_idx.i1q,:),'-', ...
         t,x2(state_idx.i1d,:),'--', ...
         t,x2(state_idx.i1q,:),'--','LineWidth',1)
    legend('$i_{1d}$','$i_{1q}$','$i_{1d}^{(2)}$','$i_{1q}^{(2)}$')
    grid on; title('Converter Currents'); axis padded

    subplot(3,1,2)
    plot(t,x1(state_idx.i2d,:),'-', ...
         t,x1(state_idx.i2q,:),'-', ...
         t,x2(state_idx.i2d,:),'--', ...
         t,x2(state_idx.i2q,:),'--','LineWidth',1)
    legend('$i_{2d}$','$i_{2q}$','$i_{2d}^{(2)}$','$i_{2q}^{(2)}$')
    grid on; title('Grid Side Currents'); axis padded

    subplot(3,1,3)
    plot(t,log1.V_PCC(1,:),'-', ...
         t,log1.V_PCC(2,:),'-', ...
         t,log2.V_PCC(1,:),'--', ...
         t,log2.V_PCC(2,:),'--','LineWidth',1)
    legend('$v_d$','$v_q$','$v_d^{(2)}$','$v_q^{(2)}$')
    grid on; title('PCC Voltage'); axis padded
end

% ============================================================
% Fig 2: Active and Reactive Power
% ============================================================
if update_plots(2)
    figure(2); clf; set(gcf,'WindowStyle','docked');

    subplot(2,1,1)
    plot(t,log1.P_conv,'-', t,log2.P_conv,'--', ...
         t,log1.P_grid,'-', t,log2.P_grid,'--', ...
         t,log1.P_PCC,'-', t,log2.P_PCC,'--','LineWidth',1)
    legend('$P_{conv}$','$P_{conv}^{(2)}$', ...
           '$P_{grid}$','$P_{grid}^{(2)}$', ...
           '$P_{PCC}$','$P_{PCC}^{(2)}$')
    grid on; title('Active Power'); axis padded

    subplot(2,1,2)
    plot(t,log1.Q_conv,'-', t,log2.Q_conv,'--', ...
         t,log1.Q_grid,'-', t,log2.Q_grid,'--', ...
         t,log1.Q_PCC,'-', t,log2.Q_PCC,'--','LineWidth',1)
    legend('$Q_{conv}$','$Q_{conv}^{(2)}$', ...
           '$Q_{grid}$','$Q_{grid}^{(2)}$', ...
           '$Q_{PCC}$','$Q_{PCC}^{(2)}$')
    grid on; title('Reactive Power'); axis padded
end

% ============================================================
% Fig 3: Synchronization States
% ============================================================
if update_plots(3)
    figure(3); clf; set(gcf,'WindowStyle','docked');

    subplot(2,1,1)
    plot(t,log1.delta_g_shifted,'-', ...
         t,log2.delta_g_shifted,'--', ...
         t,log1.delta_conv,'-', ...
         t,log2.delta_conv,'--','LineWidth',1.2)
    legend('$\delta_g$','$\delta_g^{(2)}$', ...
           '$\delta_{conv}$','$\delta_{conv}^{(2)}$')
    grid on; axis padded

    subplot(2,1,2)
    plot(t,log1.omega_g,'-', ...
         t,log2.omega_g,'--', ...
         t,log1.omega_conv,'-', ...
         t,log2.omega_conv,'--','LineWidth',1.2)
    legend('$\omega_g$','$\omega_g^{(2)}$', ...
           '$\omega_{conv}$','$\omega_{conv}^{(2)}$')
    grid on; axis padded
    title('Synchronization States'); xlabel('Time (s)')
end

% ============================================================
% Fig 4: Vref and PCC Voltage
% ============================================================
if update_plots(4)
    figure(4); clf; set(gcf,'WindowStyle','docked');

    plot(t,log1.Vrefd,'-', ...
         t,log2.Vrefd,'--', ...
         t,log1.V_PCC_mag,'-', ...
         t,log2.V_PCC_mag,'--','LineWidth',1.2)
    grid on
    title('Voltage Reference and PCC Voltage')
    xlabel('Time (s)'); ylabel('Voltage (pu)')
    legend('$V_{ref,d}$','$V_{ref,d}^{(2)}$', ...
           '$V_{PCC}$','$V_{PCC}^{(2)}$')
    axis padded
end

% ============================================================
% Fig 5: Reference Current Magnitude
% ============================================================
if update_plots(5)
    figure(5); clf; set(gcf,'WindowStyle','docked');

    plot(t,log1.Iref_mag_lim,'-', ...
         t,log2.Iref_mag_lim,'--','LineWidth',1.2)
    grid on
    title('Reference Current Magnitude')
    xlabel('Time (s)'); ylabel('$I_{ref}$ (pu)')
    legend('$I_{ref}$','$I_{ref}^{(2)}$')
    axis padded
end

% ============================================================
% Fig 6: Current Magnitudes
% ============================================================
if update_plots(6)
    figure(6); clf; set(gcf,'WindowStyle','docked');

    plot(t,log1.i1_mag,'-', t,log2.i1_mag,'--', ...
         t,log1.i2_mag,'-', t,log2.i2_mag,'--','LineWidth',1.2)
    grid on
    title('Current Magnitudes')
    xlabel('Time (s)'); ylabel('Current (pu)')
    legend('$|i_1|$','$|i_1|^{(2)}$', ...
           '$|i_2|$','$|i_2|^{(2)}$')
    axis padded
end

% ============================================================
% Fig 7: Converter Internal Voltage
% ============================================================
if update_plots(7)
    figure(7); clf; set(gcf,'WindowStyle','docked');

    subplot(2,1,1)
    plot(t,log1.Econv(1,:),'-', ...
         t,log1.Econv(2,:),'-', ...
         t,log2.Econv(1,:),'--', ...
         t,log2.Econv(2,:),'--','LineWidth',1.2)
    legend('$E_{conv,d}$','$E_{conv,q}$', ...
           '$E_{conv,d}^{(2)}$','$E_{conv,q}^{(2)}$')
    grid on; title('Converter Internal Voltage (d/q)')
    xlabel('Time (s)'); axis padded

    subplot(2,1,2)
    plot(t,log1.Econv_mag,'-', ...
         t,log2.Econv_mag,'--','LineWidth',1.2)
    grid on
    title('Converter Internal Voltage Magnitude')
    xlabel('Time (s)'); ylabel('Voltage (pu)')
    legend('$|E_{conv}|$','$|E_{conv}|^{(2)}$')
    axis padded
end

% ============================================================
% Optional: Error norm (VERY useful)
% ============================================================
figure(99); clf; set(gcf,'WindowStyle','docked');

err = x1 - x2;
err_norm = vecnorm(err,2,1);

plot(t,err_norm,'LineWidth',1.5)
grid on
title('State Error Norm')
xlabel('Time (s)')
ylabel('||x_1 - x_2||_2')



