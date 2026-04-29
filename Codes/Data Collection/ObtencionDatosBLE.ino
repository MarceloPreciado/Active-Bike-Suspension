#include <ArduinoBLE.h>
#include <Arduino_BMI270_BMM150.h>

#define SCIENCE_KIT_UUID(val) ("555a0002-" val "-467a-9538-01f0652c74e8")
#define DEBUG 0

BLEService service(SCIENCE_KIT_UUID("0000"));
BLECharacteristic accelerationCharacteristic(SCIENCE_KIT_UUID("0011"), BLENotify, 3 * sizeof(float));

String name;
unsigned long lastUpdate = 0;
const unsigned long intervalMs = 50;   // 20 Hz

void setup() {
#if DEBUG
  Serial.begin(115200);
  Serial.println("Iniciando...");
#endif

  if (!IMU.begin()) {
#if DEBUG
    Serial.println("Error al iniciar IMU");
#endif
    while (1);
  }

  if (!BLE.begin()) {
#if DEBUG
    Serial.println("Error al iniciar BLE");
#endif
    while (1);
  }

  String address = BLE.address();
  address.toUpperCase();

  name = "BLE Sense - ";
  name += address[address.length() - 5];
  name += address[address.length() - 4];
  name += address[address.length() - 2];
  name += address[address.length() - 1];

  BLE.setLocalName(name.c_str());
  BLE.setDeviceName(name.c_str());
  BLE.setAdvertisedService(service);

  service.addCharacteristic(accelerationCharacteristic);
  BLE.addService(service);

  BLE.advertise();

#if DEBUG
  Serial.print("Dispositivo BLE: ");
  Serial.println(name);
#endif
}

void loop() {
  BLE.poll();

  if (BLE.connected()) {
    if (millis() - lastUpdate >= intervalMs) {
      lastUpdate = millis();

      if (accelerationCharacteristic.subscribed() && IMU.accelerationAvailable()) {
        float x, y, z;
        float acceleration[3];

        IMU.readAcceleration(x, y, z);

        acceleration[0] = x;
        acceleration[1] = y;
        acceleration[2] = z;

        accelerationCharacteristic.writeValue((byte*)acceleration, sizeof(acceleration));
      }
    }
  }
}
