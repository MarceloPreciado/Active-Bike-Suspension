#include "Arduino_BMI270_BMM150.h"
#include <Wire.h>
#include <SPI.h>
#include <SD.h>
#include <math.h> 

// ===================================================
// CONFIGURACIÓN DE TELEMETRÍA
// ===================================================
const int PIN_CS_SD = 10;
File archivoDatos;
char nombreArchivo[16]; 

char ramBuffer[8000]; 
int ramBufferPos   = 0;
int lineasEnBuffer = 0;
const int BUFFER_LINES = 200;

int estado = 0;
unsigned long relojEstado         = 0;
unsigned long ultimoLatido        = 0;
unsigned long ultimoFlush         = 0;
unsigned long intervalMicros      = 1000; // FIJO A 1000 Hz
unsigned long ultimaMuestraUs     = 0;
unsigned long tiempoInicioGrabacion = 0;

// Analíticas en vivo
unsigned long muestrasPorSegundo = 0;
float maxGForceSegundo = 0.0;

// ───────────────────────────────────────────
void configurarBMI270() {
  Wire1.beginTransmission(0x68); Wire1.write(0x7C); Wire1.write(0x00); Wire1.endTransmission();
  delay(20);
  Wire1.beginTransmission(0x68); Wire1.write(0x40); Wire1.write(0xAC); Wire1.endTransmission();
  Wire1.beginTransmission(0x68); Wire1.write(0x41); Wire1.write(0x03); Wire1.endTransmission();
  delay(50);
}

void volcarBufferSD(bool forzarFlush) {
  if (lineasEnBuffer == 0) return;
  archivoDatos.write(ramBuffer, ramBufferPos);
  ramBufferPos   = 0;
  lineasEnBuffer = 0;
  if (forzarFlush) archivoDatos.flush();
}

// ───────────────────────────────────────────
void setup() {
  Serial.begin(115200);
  
  // Arranque inteligente: Espera 3 seg a la PC. Si no hay PC, arranca con batería.
  unsigned long tInicio = millis();
  while (!Serial && millis() - tInicio < 3000);

  Serial.println("=========================================");
  Serial.println("🚴 TELEMETRÍA V4 - CORTES DE 10 SEGUNDOS 🚴");
  Serial.println("=========================================");

  if (!SD.begin(PIN_CS_SD)) {
    if(Serial) Serial.println("❌ ERROR FATAL: SD falló."); estado = 99; return;
  }

  // Generador automático de archivos
  for (int i = 0; i < 1000; i++) {
    snprintf(nombreArchivo, sizeof(nombreArchivo), "run_%03d.csv", i);
    if (!SD.exists(nombreArchivo)) break;
  }
  
  archivoDatos = SD.open(nombreArchivo, FILE_WRITE);
  if (!archivoDatos) {
    if(Serial) Serial.println("❌ ERROR FATAL: No se pudo crear archivo."); estado = 99; return;
  }
  
  if(Serial) { Serial.print("✅ Archivo creado: "); Serial.println(nombreArchivo); }
  
  archivoDatos.println("timestamp,x,y,z");
  archivoDatos.flush();

  if (!IMU.begin()) {
    if(Serial) Serial.println("❌ ERROR FATAL: IMU falló."); estado = 99; return;
  }
  
  Wire1.setClock(400000);
  configurarBMI270();
  if(Serial) Serial.println("✅ IMU a 1000 Hz / 16g.");

  estado = 1;
  relojEstado = millis();
  if(Serial) Serial.println("\n>>> 🟢 TIENES 10 SEGUNDOS DE PREPARACIÓN <<<");
}

// ───────────────────────────────────────────
void loop() {
  unsigned long ahoraMs = millis();
  unsigned long ahoraUs = micros();

  // ── ESTADO 1: Cuenta regresiva (10 seg) ──
  if (estado == 1) {
    if (ahoraMs - ultimoLatido >= 1000) {
      int segs = 10 - (int)((ahoraMs - relojEstado) / 1000);
      if(Serial) { Serial.print("⏳ Inicio en "); Serial.print(segs); Serial.println(" seg..."); }
      ultimoLatido = ahoraMs;
    }

    if (ahoraMs - relojEstado >= 10000) {
      tiempoInicioGrabacion = ahoraMs;
      ultimaMuestraUs       = ahoraUs;
      ultimoLatido          = ahoraMs;
      ultimoFlush           = ahoraMs;
      muestrasPorSegundo    = 0;
      maxGForceSegundo      = 0.0;
      
      if(Serial) Serial.println("\n🔥 >>> GRABANDO (Límite: 10 Segundos) <<< 🔥");
      estado = 2;
    }
    return;
  }

  // ── ESTADO 2: Grabación Activa ──
  if (estado == 2) {
    
    unsigned long tiempoGrabado = ahoraMs - tiempoInicioGrabacion;

    // 🛑 EL CORTADOR AUTOMÁTICO DE 10 SEGUNDOS 🛑
    if (tiempoGrabado >= 10000) {
      volcarBufferSD(true); // Vacía cualquier dato que haya quedado en la memoria RAM
      archivoDatos.close(); // Cierra el archivo de forma 100% segura
      
      if(Serial) {
        Serial.println("\n=========================================");
        Serial.println("🏁 GRABACIÓN DE 10 SEGUNDOS COMPLETADA 🏁");
        Serial.print("✅ Archivo guardado y cerrado: "); Serial.println(nombreArchivo);
        Serial.println("👉 Presiona el botón RESET en la placa para grabar otro archivo.");
        Serial.println("=========================================");
      }
      
      estado = 3; // Pasamos a estado "Dormido"
      return;
    }

    // Resumen de analíticas cada segundo
    if (ahoraMs - ultimoLatido >= 1000) {
      if(Serial) {
        Serial.print("⏱️ Seg: "); Serial.print(tiempoGrabado / 1000);
        Serial.print(" | 📈 Muestras: "); Serial.print(muestrasPorSegundo);
        Serial.print(" | 💥 Impacto Máx: "); Serial.print(maxGForceSegundo); Serial.println(" G");
      }
      muestrasPorSegundo = 0;
      maxGForceSegundo = 0.0;
      ultimoLatido = ahoraMs;
    }

    if (ahoraMs - ultimoFlush >= 2000) {
      volcarBufferSD(true);
      ultimoFlush = ahoraMs;
    }

    if (ahoraUs - ultimaMuestraUs < intervalMicros) return;
    ultimaMuestraUs += intervalMicros;

    Wire1.beginTransmission(0x68); Wire1.write(0x0C); Wire1.endTransmission(false);
    Wire1.requestFrom((uint8_t)0x68, (uint8_t)6);

    if (Wire1.available() >= 6) {
      int16_t x_raw = (int16_t)(Wire1.read() | (Wire1.read() << 8));
      int16_t y_raw = (int16_t)(Wire1.read() | (Wire1.read() << 8));
      int16_t z_raw = (int16_t)(Wire1.read() | (Wire1.read() << 8));

      float x_ms2 = (x_raw / 2048.0) * 9.81;
      float y_ms2 = (y_raw / 2048.0) * 9.81;
      float z_ms2 = (z_raw / 2048.0) * 9.81;

      float g_actual = sqrt((x_raw/2048.0)*(x_raw/2048.0) + (y_raw/2048.0)*(y_raw/2048.0) + (z_raw/2048.0)*(z_raw/2048.0));
      if (g_actual > maxGForceSegundo) maxGForceSegundo = g_actual;

      muestrasPorSegundo++;

      int n = snprintf(ramBuffer + ramBufferPos, sizeof(ramBuffer) - ramBufferPos, "%lu,%.2f,%.2f,%.2f\n", tiempoGrabado, x_ms2, y_ms2, z_ms2);
      if (n > 0) { ramBufferPos += n; lineasEnBuffer++; }
      if (lineasEnBuffer >= BUFFER_LINES) volcarBufferSD(false);
    }
    return;
  }

  // ── ESTADO 3: Grabación Terminada (Dormido) ──
  if (estado == 3) {
    // La placa ya no hace nada. Espera a que la reinicies o la desconectes.
    delay(1000);
  }

  // ── ESTADO 99: Error fatal ──
  if (estado == 99) {
    delay(100);
  }
}