#!/usr/bin/env python3

import zmq
import struct
from datetime import datetime


__all__ = ["VRTSubscriber"]


class VRTSubscriber:
    def __init__(self, host="127.0.0.1", port=50100, hwm=10000, timeout=1):
        context = zmq.Context()
        self.zmq_subscriber = context.socket(zmq.SUB)
        self.zmq_subscriber.setsockopt(zmq.RCVHWM, hwm)
        self.zmq_subscriber.setsockopt_string(zmq.SUBSCRIBE, "")
        self.zmq_subscriber.setsockopt(zmq.RCVTIMEO, timeout * 1000)
        connect_string = f"tcp://{host}:{port}"
        self.zmq_subscriber.connect(connect_string)
        self.timeout = timeout

    def get_context_packet(self):
        """Get the first context package from the stream (discard data packets)"""
        while True:
            try:
                packet = self.zmq_subscriber.recv()
            except zmq.Again as e:
                raise RuntimeError("Timed out connecting to VRT stream") from e
            if packet[0] == 73:
                return packet

    def get_time(self):
        """
        Get the timestamp from the stream (currently only with integer seconds)
        """
        context_packet = self.get_context_packet()
        (unixtime_seconds,) = struct.unpack("!I", context_packet[16:20])
        return datetime.utcfromtimestamp(unixtime_seconds)

    def is_live(self, tolerance_s=2):
        """
        Returns True if the VRT stream is running now; False if the timestamp
        in the stream does not match the current system clock.
        """
        stream_time = self.get_time()
        now = datetime.now()
        return abs((now - stream_time).seconds) < tolerance_s + 0.1

    def get_frequency(self):
        """Get the center frequency from the stream"""
        raise NotImplementedError()

    def get_bandwidth(self):
        """Get the bandwidth from the stream"""
        raise NotImplementedError()


if __name__ == "__main__":
    client = VRTSubscriber()
    print(client.get_time())
