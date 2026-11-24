"""Transport-stub checks for client safety; no EOS source-data fixtures."""

import unittest
import os
from types import SimpleNamespace

import numpy as np

from api_client import EOSClient, validate_success


class ClientTests(unittest.TestCase):
    def setUp(self):
        self.result = np.ones(125)
        self.result[:4] = [1, 0, 1, 1]
        self.result[23:31] = 0
        self.result[25] = 1
        self.client = object.__new__(EOSClient)
        self.client.query = lambda *arguments: self.result.copy()

    def test_valid_state(self):
        state = self.client.state("candidate", "component", 0., 0.)
        self.assertTrue(all(np.isfinite(value) for value in state.values()))

    def test_missing_candidate_callback_is_rejected(self):
        self.result[2] = 0
        with self.assertRaises(ValueError):
            self.client.state("candidate", "component", 0., 0.)

    def test_unexpected_legacy_callback_is_rejected(self):
        with self.assertRaises(ValueError):
            self.client.state("legacy", "component", 0., 0.)

    def test_invalid_fractions_are_rejected(self):
        for fraction in (-.1, .5, 1.1, np.nan):
            with self.subTest(fraction=fraction):
                self.result[25] = fraction
                with self.assertRaises(ValueError):
                    self.client.state("candidate", "component", 0., 0.)

    def test_wrong_component_identity_is_rejected(self):
        self.result[25] = 0
        self.result[23] = 1
        self.result[2] = 0
        with self.assertRaises(ValueError):
            self.client.state("candidate", "component", 0., 0.)
        self.assertTrue(np.isfinite(self.client.state("candidate", "final", 0., 0.)["P"]))

    def test_reconstructed_overflow_is_rejected(self):
        self.result[4] = 1000
        with self.assertRaises(ValueError):
            self.client.state("candidate", "component", 0., 0.)

    def test_hybrid_failure_is_rejected(self):
        self.result[3] = 0
        with self.assertRaises(ValueError):
            self.client.state("candidate", "component", 0., 0.)

    def test_protocol_shape_and_eos_errors_are_rejected(self):
        with self.assertRaises(ValueError):
            validate_success("candidate", "component", self.result[:-1])
        self.result[1] = -7999
        with self.assertRaises(ValueError):
            self.client.state("candidate", "component", 0., 0.)

    def pipe_response(self, payload):
        read_descriptor, write_descriptor = os.pipe()
        self.client.timeout = .02
        try:
            if payload is not None:
                os.write(write_descriptor, payload)
            if payload == b"":
                os.close(write_descriptor)
                write_descriptor = None
            with os.fdopen(read_descriptor, "r") as stream:
                self.client.process = SimpleNamespace(stdout=stream)
                return self.client.read_response()
        finally:
            if write_descriptor is not None:
                os.close(write_descriptor)

    def test_response_deadline(self):
        with self.assertRaises(TimeoutError):
            self.pipe_response(None)

    def test_response_eof(self):
        with self.assertRaises(RuntimeError):
            self.pipe_response(b"")

    def test_response_line(self):
        self.assertEqual(self.pipe_response(b"1 0 1 1\n"), "1 0 1 1\n")


if __name__ == "__main__":
    unittest.main()
