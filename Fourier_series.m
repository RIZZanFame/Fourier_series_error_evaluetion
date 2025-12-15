% Tale programma, data una funzione nota, mi calcola, sino all'N armonica, 
% la sua serie di Fourier e ne valuta l'errore.

% Funzione y(x) in -pi<x<pi

% definisco il quantitativo di campioni e la funzione in esame
x = linspace(-pi,pi,200000);
y = abs(sin(x));

% definisco il numero delle armoniche e calcolo i 3 coeficenti
n = (1:1:200)';
a_0 = (1/(2*pi)) .* trapz(x, y);
a_k = (1/pi) .* trapz(x, y .* cos(n.*x), 2);
b_k = (1/pi) .* trapz(x, y .* sin(n.*x), 2);

% computo la serie di Fourier con la sua formula canonica
fs = a_0 + (a_k.' * cos(n.*x)) + (b_k.' * sin(n.*x));

% creo il vettore delle 'distanze'
norma_2_coll = [];
coincidence_x_coll = [];
k = 1;
toll = sqrt(eps);

for ii = 1:length(fs)
    norma_2_coll(ii) = abs(y(ii) - fs(ii));
    if norma_2_coll(ii) < toll
        coincidence_x_coll(k) = ii;
        k = k + 1;
    end
end

% plotto il tutto (tanto per anche questo)
figure
plot(x, y, 'LineWidth', 1.5); hold on
plot(x, fs, '--', 'LineWidth', 1.5);

grid on
legend('y(x)', 'Serie di Fourier')
xlabel('x')
ylabel('valore')

figure
plot(x, norma_2_coll, 'LineWidth', 1.5); hold on

plot(x(coincidence_x_coll), ...
     norma_2_coll(coincidence_x_coll), ...
     'ro', 'MarkerSize', 5, 'LineWidth', 1.2)

grid on
legend('Funzione degli errori', 'Coincidenze numeriche')
xlabel('x')
ylabel('Errore')

