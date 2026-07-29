"""Common natural-language router for specialized py4siesta agents."""

from py4siesta_agent.shared.types import RouteDecision


class AgentRouter:
    """Route a request using only declared specialized-agent capabilities."""

    def __init__(self, model_client, agents):
        self.model_client = model_client
        self.agents = dict(agents)

    def route(self, request):
        capabilities = "\n".join(
            "- %s: %s" % (name, agent.capability)
            for name, agent in sorted(self.agents.items())
        )
        decision = self.model_client.structured(
            (
                "Select an available specialized agent only when the request matches "
                "its declared capability. Otherwise select unsupported. Available:\n"
                + capabilities
            ),
            request,
            RouteDecision,
        )
        if decision.agent == "unsupported" or decision.agent not in self.agents:
            return {
                "ok": False,
                "status": "unsupported",
                "selected_agent": None,
                "reason": decision.reason,
            }
        result = self.agents[decision.agent].run(request)
        result.setdefault("selected_agent", self.agents[decision.agent].name)
        result.setdefault("routing_decision", decision.model_dump())
        return result
