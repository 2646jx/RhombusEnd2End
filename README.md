# Rhombus: Fast Homomorphic Matrix-Vector Multiplication for Secure Two-Party Inference

This repository is developed based on [OpenCheetah](https://github.com/Alibaba-Gemini-Lab/OpenCheetah), and provides a proof-of-concept implementation of the paper [Rhombus](https://eprint.iacr.org/2024/1611.pdf). In summary, this project contains the following components:

- Tests for the protocols: matrix-vector multiplication, matrix-matrix multiplication
- End-to-end implementation of secure two-party inference of ResNet50 model

If someone wants to learn more details about the algorithms in [Rhombus](https://eprint.iacr.org/2024/1611.pdf), refer to [here](https://github.com/2646jx/Rhombus), in which we provide a pure implementation without the network transmission.

## Building

Run the following command to build the dependencies:

```PowerShell
bash scripts/build-deps.sh
```

After building the dependencies successfully, run the following command to build this project:

```PowerShell
bash scripts/build.sh
```

This will generate four executable files.

## Testing

Before performing the tests, execute the script file `scripts/throttle.sh` to mimic the network environment. For example,

Run the following command to mimic an LAN setting (3Gbps, 0.3ms):

```PowerShell
bash scripts/throttle.sh lan
```

Run the following command to mimic a WAN setting (100Mbps, 40ms):

```PowerShell
bash scripts/throttle.sh wan
```

Run the following command to remove the network traffic restriction:

```PowerShell
bash scripts/throttle.sh del
```

After setting the network environment, launch two terminals, representing the client and server, respectively.

### Module test

To test matrix-vector multiplication protocol, execute the following command in the client terminal:

```PowerShell
./build/bin/rhombus_matvec 2 12345
```

Correspondingly, execute the following command in the server terminal:

```PowerShell
./build/bin/rhombus_matvec 1 12345
```

The matrix-matrix multiplication protocol can be tested in a similar way as above.

In client terminal:

```PowerShell
./build/bin/rhombus_matmul 2 12345
```

In server terminal:

```PowerShell
./build/bin/rhombus_matmul 1 12345
```

After running the protocols successfully, the performance data will be printed. The total communication volume of the protocol is the sum of the client's and server's transmitted data. To configure the parameters of the protocol, e.g., the dimensions of the matrices, the number of threads, you can modify the parameters in the  corresponding files, then recompile them.

### End-to-end test for ResNet50

To run [Cheetah](https://eprint.iacr.org/2022/207.pdf), execute the following command in client terminal:

```PowerShell
sh scripts/run-client.sh cheetah resnet50
```

In the server terminal, execute:

```PowerShell
sh scripts/run-server.sh cheetah resnet50
```

To run [Rhombus](https://eprint.iacr.org/2024/1611.pdf), you can directly change the `cheetah` by `rhombus` in the above commands. In particular, run this
command in client terminal:

```PowerShell
sh scripts/run-client.sh rhombus resnet50
```

In server terminal, run:

```PowerShell
sh scripts/run-server.sh cheetah resnet50
```

After running the two programs, the performance results will be printed in the server's terminal.

## LICENSE

This project is licensed under the MIT License. See the [LICENSE] file for details.