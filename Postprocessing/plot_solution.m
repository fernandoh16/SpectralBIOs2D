clear all
close all
clc

M = load('solution.txt');
figure
plot(M(:,2))
hold on
plot(M(:,3))
dh = M(2,1)-M(1,1);
Error = norm(M(:,2)-M(:,3))*sqrt(dh);