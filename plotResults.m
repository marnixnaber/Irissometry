
% plot basic irissometry results
h = figure();
subplot(4,1,1);
plot(eyeData(:,6),eyeData(:,2))
xlabel('Time [s]')
ylabel([eyeDataHeader{2} ' [pixels]'])
title(videoFileNames{vidNum})

subplot(4,1,2);
plot(eyeData(:,6),eyeData(:,3),'r')
hold on
plot(eyeData(:,6),eyeData(:,30),'m')
plot(eyeData(:,6),eyeData(:,31),'r:')
plot(eyeData(:,6),eyeData(:,32),'m:')
xlabel('Time [s]')
ylabel(['Coordinate [pixels]'])
% legend('Pupil','Iris','Inner iris','Outer iris')
legend({eyeDataHeader{[3 30:32]}})

subplot(4,1,3);
plot(eyeData(:,6),eyeData(:,4),'b')
hold on
plot(eyeData(:,6),eyeData(:,43),'c')
plot(eyeData(:,6),eyeData(:,44),'b:')
plot(eyeData(:,6),eyeData(:,45),'c:')
xlabel('Time [s]')
ylabel('Y coordinate [pixels]')
legend({eyeDataHeader{[4 43:45]}})

subplot(4,1,4);
plot(eyeData(:,6),eyeData(:,5))
xlabel('Time [s]')
ylabel(eyeDataHeader{5})


% nanstd(eyeData(:,[3 30:32]))
% nanstd(eyeData(:,[4 43:45]))