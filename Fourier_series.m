% Tale programma, data una funzione nota, mi calcola, sino all'N armonica, 
% la sua serie di Fourier e ne valuta l'errore.

% Funzione y(x) in -pi<x<pi

% definisco il quantitativo di campioni e la funzione in esame
x = linspace(-pi,pi,200000);
y = log(abs(x)) + abs(sin(x)) + abs(cos(x));

% definisco il numero delle armoniche e calcolo i 3 coeficenti
n = (1:1:1000)';
a_0 = (1/(2*pi)) .* trapz(x, y);
a_k = (1/pi) .* trapz(x, y .* cos(n.*x), 2);
b_k = (1/pi) .* trapz(x, y .* sin(n.*x), 2);

% computo la serie di Fourier con la sua formula canonica
fs = a_0 + (a_k.' * cos(n.*x)) + (b_k.' * sin(n.*x));

% calcolo l'errore (tanto per)
err = y - fs;
norm_L2 = sqrt(trapz(x, err.^2));

% plotto il tutto (tanto per anche questo)
figure
plot(x, y, 'LineWidth', 1.5); hold on
plot(x, fs, '--', 'LineWidth', 1.5);
grid on
legend('y(x)', 'Serie di Fourier')
xlabel('x')
ylabel('valore')