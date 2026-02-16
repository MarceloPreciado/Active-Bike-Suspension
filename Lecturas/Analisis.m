%% Accelerometer FFT Analyzer (robusto con celdas vacías)
clear; close all; clc;

cfg.windowType      = "hann";   % "hann","hamming","blackman","rect"
cfg.useZeroPad      = true;
cfg.zeroPadFactor   = 4;
cfg.maxFreqPlotHz   = [];       % [] = hasta Nyquist
cfg.plotTimeSignals = true;
cfg.plotFFT         = true;

%% Seleccionar archivos
[files, path] = uigetfile({'*.csv','CSV Files (*.csv)'}, ...
    'Selecciona uno o varios CSV', 'MultiSelect','on');
if isequal(files,0)
    disp("No seleccionaste archivos. Saliendo...");
    return;
end
if ischar(files); files = {files}; end

for k = 1:numel(files)
    filename = fullfile(path, files{k});
    fprintf("\n=== Procesando: %s ===\n", files{k});

    % --- Lectura simple (sin WhitespaceRule) ---
    opts = detectImportOptions(filename);
    T = readtable(filename, opts);

    % Normalizar nombres
    varNames = string(T.Properties.VariableNames);
    varNorm  = normalizeVarNames(varNames);

    % Buscar columnas (tu caso: timestamp, LinearAccc, AccX, AccY, AccZ)
    idx.t   = findFirstMatch(varNorm, ["timestamp","time","tiempo","t"]);
    idx.lin = findFirstMatch(varNorm, ["linearaccc","linearacc","linearaccel","aceleracionlineal","linealacc","linear"]);
    idx.ax  = findFirstMatch(varNorm, ["accx","accelx","ax"]);
    idx.ay  = findFirstMatch(varNorm, ["accy","accely","ay"]);
    idx.az  = findFirstMatch(varNorm, ["accz","accelz","az"]);

    % Si algo falta, pedir índice
    fields = ["t","lin","ax","ay","az"];
    labels = ["Timestamp/tiempo","Linear Acc","AccX","AccY","AccZ"];
    for i = 1:numel(fields)
        if isempty(idx.(fields(i)))
            disp("Columnas disponibles:");
            for j = 1:numel(varNames)
                fprintf("  %2d) %s\n", j, varNames(j));
            end
            idx.(fields(i)) = input("Indice para " + labels(i) + ": ");
        end
    end

    % Extraer (permitiendo NaN en señales)
    tRaw = T{:, idx.t};
    lin  = T{:, idx.lin};
    ax   = T{:, idx.ax};
    ay   = T{:, idx.ay};
    az   = T{:, idx.az};

    % Asegurar columna
    tRaw = tRaw(:); lin = lin(:); ax = ax(:); ay = ay(:); az = az(:);

    % Quitar filas donde el tiempo sea NaN
    keepT = isfinite(tRaw);
    tRaw = tRaw(keepT);
    lin  = lin(keepT);
    ax   = ax(keepT);
    ay   = ay(keepT);
    az   = az(keepT);

    % Ordenar por tiempo
    [tRaw, sidx] = sort(tRaw);
    lin = lin(sidx); ax = ax(sidx); ay = ay(sidx); az = az(sidx);

    % --- Heurística de timestamp (epoch ms vs epoch s vs ya en s) ---
    % epoch en ms ~ 1.7e12, epoch en s ~ 1.7e9
    if mean(tRaw, "omitnan") > 1e11
        tSec = tRaw / 1000;  % ms -> s
        disp("Detecté timestamp tipo epoch en ms. Convirtiendo a segundos...");
    else
        tSec = tRaw;         % ya está en s
    end

    % Hacer tiempo relativo desde 0
    tSec = tSec - tSec(1);

    % Procesar por señal (porque hay muchos huecos)
    signals = {lin, ax, ay, az};
    sigNames = ["LinearAcc","AccX","AccY","AccZ"];

    if cfg.plotTimeSignals
        figure('Name', "Time Signals - " + files{k}, 'Color','w');
        tiledlayout(4,1,'Padding','compact','TileSpacing','compact');
    end

    if cfg.plotFFT
        figure('Name', "FFT - " + files{k}, 'Color','w');
        tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
    end

    for i = 1:4
        sig = signals{i};

        % Usar SOLO filas donde esta señal exista (no NaN)
        m = isfinite(tSec) & isfinite(sig);
        ti = tSec(m);
        xi = sig(m);

        if numel(xi) < 8
            warning("%s: Muy pocos puntos válidos para FFT (%d).", sigNames(i), numel(xi));
            continue;
        end

        % Estimar Fs desde esta señal (no desde el CSV completo)
        dt = median(diff(ti));
        Fs = 1/dt;

        % Preprocesado básico
        xi = xi - mean(xi, 'omitnan');
        xi = detrend(xi);

        % Plot tiempo
        if cfg.plotTimeSignals
            figure(findobj('Name',"Time Signals - " + files{k}));
            nexttile(i);
            plot(ti, xi); grid on;
            ylabel(sigNames(i));
            if i==1
                title("Señales en el tiempo (solo puntos válidos) - " + files{k}, 'Interpreter','none');
            end
            if i==4
                xlabel("Tiempo [s]");
            end
        end

        % FFT
        [f, A] = singleSidedFFT(xi, Fs, cfg.windowType, cfg.useZeroPad, cfg.zeroPadFactor);

        if cfg.plotFFT
            figure(findobj('Name',"FFT - " + files{k}));
            nexttile(i);
            plot(f, A); grid on;
            title(sigNames(i) + " | Fs≈" + num2str(Fs, '%.2f') + " Hz");
            xlabel("Frecuencia [Hz]");
            ylabel("Amplitud");
            if ~isempty(cfg.maxFreqPlotHz)
                xlim([0 cfg.maxFreqPlotHz]);
            else
                xlim([0 Fs/2]);
            end
        end
    end
end

%% ===== Funciones =====
function out = normalizeVarNames(names)
    out = lower(names);
    out = replace(out, [" ", "-", "_", "(", ")", "[", "]", "{", "}", ".", ","], "");
    out = replace(out, ["á","é","í","ó","ú","ñ"], ["a","e","i","o","u","n"]);
end

function idx = findFirstMatch(varNorm, candidates)
    idx = [];
    for c = candidates
        hits = find(varNorm == c, 1, 'first');
        if ~isempty(hits)
            idx = hits;
            return;
        end
    end
end

function [f, A] = singleSidedFFT(x, Fs, windowType, useZeroPad, zeroPadFactor)
    x = x(:);
    N = length(x);

    w = getWindow(windowType, N);
    xw = x .* w;

    cg = mean(w); % coherent gain

    if useZeroPad
        Nfft = 2^nextpow2(zeroPadFactor * N);
    else
        Nfft = 2^nextpow2(N);
    end

    X = fft(xw, Nfft);

    P2 = abs(X / (N * cg));
    P1 = P2(1:Nfft/2+1);
    P1(2:end-1) = 2*P1(2:end-1);

    f = Fs*(0:(Nfft/2))/Nfft;
    A = P1;
end

function w = getWindow(windowType, N)
    switch lower(string(windowType))
        case "hann"
            w = hann(N);
        case "hamming"
            w = hamming(N);
        case "blackman"
            w = blackman(N);
        case "rect"
            w = ones(N,1);
        otherwise
            w = hann(N);
    end
end
