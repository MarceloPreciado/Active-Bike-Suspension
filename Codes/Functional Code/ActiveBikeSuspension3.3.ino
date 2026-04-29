#include <SuspensionBike3.2_inferencing.h>
#include <Arduino_BMI270_BMM150.h>
#include <ArduinoBLE.h>
#include <Servo.h> 

/* BLE -------------------------------------------------------------------- */
BLEService terrainService("180C");
BLEStringCharacteristic terrainCharacteristic("2A56", BLERead | BLENotify, 50);

/* Servo ------------------------------------------------------------------ */
Servo miServo;           
#define SERVO_PIN 9    

/* Modelo ------------------------------------------------------------------ */
#define CONVERT_G_TO_MS2    9.80665f
#define MAX_ACCEPTED_RANGE  2.0f

static float features[EI_CLASSIFIER_DSP_INPUT_FRAME_SIZE];
static size_t feature_ix = 0;

/* Votación por mayoría ---------------------------------------------------- */
#define INTERVALO_ENVIO_MS  2000

static unsigned long tiempoVentana = 0;
static int votos[EI_CLASSIFIER_LABEL_COUNT] = {0};

/* Prototipos --------------------------------------------------------------- */
void enviarClasificacionBLE(const char* label, float value);
int obtenerIndiceGanador(ei_impulse_result_t *result);
int obtenerEtiquetaMasVotada();
void reiniciarVotos();
void moverServoSegunTerreno(const char* etiqueta); 

/* Setup -------------------------------------------------------------------- */
void setup() {
    Serial.begin(115200);
    unsigned long startTime = millis();
    while (!Serial && millis() - startTime < 2000);
    delay(1000);

    Serial.println("Iniciando Nano 33 BLE Sense + Servo...");

    // 4. Configuración inicial del Servo
    miServo.attach(SERVO_PIN);
    miServo.write(90); 
    Serial.println("Servo inicializado en 90 grados (Medio).");

    if (!IMU.begin()) {
        Serial.println("Error: no se pudo iniciar IMU.");
        while (1);
    }

    if (!BLE.begin()) {
        Serial.println("Error: no se pudo iniciar BLE.");
        while (1);
    }

    BLE.setLocalName("Nano33_Terreno");
    BLE.setDeviceName("Nano33_Terreno");
    BLE.setAdvertisedService(terrainService);
    terrainService.addCharacteristic(terrainCharacteristic);
    BLE.addService(terrainService);
    terrainCharacteristic.writeValue("Esperando...");
    BLE.advertise();

    tiempoVentana = millis();
}

/* Loop --------------------------------------------------------------------- */
void loop() {
    BLE.poll();

    static bool wasConnected = false;
    BLEDevice central = BLE.central();

    if (central && !wasConnected) {
        Serial.print("Conectado a: ");
        Serial.println(central.address());
        wasConnected = true;
    }
    if (!central && wasConnected) {
        Serial.println("Central desconectada.");
        wasConnected = false;
    }

    while (feature_ix < EI_CLASSIFIER_DSP_INPUT_FRAME_SIZE) {
        float x, y, z;
        if (!IMU.accelerationAvailable()) {
            delay(1);
            return;
        }
        IMU.readAcceleration(x, y, z);
        features[feature_ix++] = x * CONVERT_G_TO_MS2;
        if (feature_ix < EI_CLASSIFIER_DSP_INPUT_FRAME_SIZE) features[feature_ix++] = y * CONVERT_G_TO_MS2;
        if (feature_ix < EI_CLASSIFIER_DSP_INPUT_FRAME_SIZE) features[feature_ix++] = z * CONVERT_G_TO_MS2;
    }

    signal_t signal;
    int err = numpy::signal_from_buffer(features, EI_CLASSIFIER_DSP_INPUT_FRAME_SIZE, &signal);
    if (err != 0) {
        feature_ix = 0;
        return;
    }

    ei_impulse_result_t result = { 0 };
    EI_IMPULSE_ERROR res = run_classifier(&signal, &result, false);
    if (res != EI_IMPULSE_OK) {
        feature_ix = 0;
        return;
    }

    int ganador = obtenerIndiceGanador(&result);
    if (ganador >= 0) {
        votos[ganador]++;
    }

    // Cada 2 segundos: Actuar según votación
    if (millis() - tiempoVentana >= INTERVALO_ENVIO_MS) {
        int masVotado = obtenerEtiquetaMasVotada();

        if (masVotado >= 0) {
            const char* etiqueta = result.classification[masVotado].label;

            Serial.print(">>> Resultado ventana: ");
            Serial.println(etiqueta);

            // 5. Mover el servo según la etiqueta ganadora
            moverServoSegunTerreno(etiqueta);

            enviarClasificacionBLE(etiqueta, (float)votos[masVotado]);
        }

        reiniciarVotos();
        tiempoVentana = millis();
    }

    feature_ix = 0;
}

/* Función para mover el Servo --------------------------------------------- */
void moverServoSegunTerreno(const char* etiqueta) {

    
    if (strcmp(etiqueta, "Rough") == 0) {
        miServo.write(180);
        Serial.println("Servo -> 180° (Antihorario)");
    } 
    else if (strcmp(etiqueta, "Liso") == 0) {
        miServo.write(0);
        Serial.println("Servo -> 0°");
    } 
    else if (strcmp(etiqueta, "Medio") == 0) {
        miServo.write(90);
        Serial.println("Servo -> 90° (Posición inicial)");
    }
}

/* Resto de funciones (Sin cambios significativos) ------------------------- */
int obtenerIndiceGanador(ei_impulse_result_t *result) {
    int best_ix = -1;
    float best_val = -1.0f;
    for (size_t ix = 0; ix < EI_CLASSIFIER_LABEL_COUNT; ix++) {
        if (result->classification[ix].value > best_val) {
            best_val = result->classification[ix].value;
            best_ix = ix;
        }
    }
    return best_ix;
}

int obtenerEtiquetaMasVotada() {
    int best_ix = -1;
    int best_votos = -1;
    for (size_t ix = 0; ix < EI_CLASSIFIER_LABEL_COUNT; ix++) {
        if (votos[ix] > best_votos) {
            best_votos = votos[ix];
            best_ix = ix;
        }
    }
    return best_ix;
}

void reiniciarVotos() {
    for (size_t ix = 0; ix < EI_CLASSIFIER_LABEL_COUNT; ix++) votos[ix] = 0;
}

void enviarClasificacionBLE(const char* label, float value) {
    terrainCharacteristic.writeValue(label);
}