# Deferred HTTP transport findings

These accepted limitations apply to the shared HTTP clients. Resolve them with
focused regression tests when extending the affected behavior.

| Finding | Current behavior and reason for deferral | Follow-up |
| --- | --- | --- |
| Malformed request headers | Backend validation can raise different exception types. Valid requests use the shared transport error contract; invalid-input normalization is deferred. | Validate malformed headers consistently before a network attempt. |
| Raw-body Content-Type | Requests and aiohttp choose different defaults when callers omit Content-Type. Callers can provide an explicit header. | Choose and test a shared default before a provider relies on implicit content typing. |
| Invalid CA paths | TLS setup can raise a configuration exception directly. This preserves the underlying configuration failure. | Define a shared configuration-error contract before normalizing these exceptions. |
| Duplicated sync and async loops | Shared policy owns budgets, retries and redirects, but the backend I/O loops remain separate. Refactoring them would expand the current change. | Keep conformance tests shared; consolidate only where it simplifies both implementations. |
| CA environment variables | Native and Requests clients honor REQUESTS_CA_BUNDLE and CURL_CA_BUNDLE even with trust_env=False; aiohttp uses the documented download SSL context. Explicit ca_bundle provides consistent trust configuration. | Decide on a shared environment policy before changing the documented behavior. |

Cache entries with a `Vary` response header are deliberately not stored. Before
adopting the transport in a provider that returns `Vary: Accept-Encoding`, decide
whether to support variant-aware caching; the current restriction is documented
in [HTTP transport](docs/io/http-transport.md).

The dependency lower-bound checker imports `packaging` through the locked
environment's transitive dependencies. If the checker must run in a minimal
environment, declare that dependency explicitly or remove the import.
